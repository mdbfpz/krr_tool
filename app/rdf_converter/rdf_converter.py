from rdflib import Graph, URIRef, Literal, Namespace
from rdflib.namespace import RDF, RDFS, XSD 
import xml.etree.ElementTree as ET
import copy

# Define namespaces
BASE = Namespace("https://aware-sesar.eu/")
FIXM = Namespace("http://www.fixm.aero/flight/4.3/")
FB = Namespace("http://www.fixm.aero/base/4.3/")
HMI = Namespace("https://aware-sesar.eu/hmi/") # TODO: do i need this or can add to base?
XS = Namespace("http://www.w3.org/2001/XMLSchema/")

class RDFConverter:
    def __init__(self):
        """Initialize the RDF graph and bind namespaces."""

        self.namespaces = {
            "": BASE,
            "fx": FIXM,
            "fb": FB,
            "hmi": HMI,
            "xs": XS,
            "xsd": XSD,
            "rdfs": RDFS,
            "rdf": RDF,
        }

        # Initialize the RDF graph and bind prefixes
        self.graph = Graph()
        for prefix, ns in self.namespaces.items():
            self.graph.bind(prefix, ns, override=True)

        self.data_repository = {}
        self.track_number_callsign_map = {} # key = track number; value = callsign
        self.uuid_callsign_map = {}  # key = flight plan UUID; value = callsign
        self.last_timestamp = None  # Store the last processed timestamp

        self.clearance_keys = ["heading"]

        self.aircraft_id = ""  # Default aircraft ID
        self.id_counter = {}  # Counter for generating unique IDs for each type

        self.to_fl_const = 3.28084 / 100

    def _generate_index(self, tag, timestamp):
        """Generates a unique index for a given tag."""
        if tag not in self.id_counter:
            self.id_counter[tag] = 1
        else:
            self.id_counter[tag] += 1
        return f"{tag}_{self.aircraft_id}_{timestamp}_{self.id_counter[tag]}"

    def _add_literal(self, subject, predicate, value, datatype=XSD.string):
        """Add a literal triple to the graph."""
        self.graph.add((subject, predicate, Literal(value, datatype=datatype)))

    def _process_position(self, parent_uri, position_data, fx_namespace, fb_namespace):
        """
        Create RDF triples for a structured position object.
        - Keys with nested dicts are processed under 'pos'
        - All other top-level keys are added under 'position'
        """
        position_uri = URIRef(f"{parent_uri}_position")
        pos_uri = URIRef(f"{position_uri}_pos")

        # Structure: parent -> position -> pos
        self.graph.add((parent_uri, fx_namespace.position, position_uri))
        self.graph.add((position_uri, RDF.type, fx_namespace.Position))
        self.graph.add((position_uri, fb_namespace.pos, pos_uri))
        self.graph.add((pos_uri, RDF.type, fb_namespace.Pos))

        for key, value in position_data.items():
            if isinstance(value, dict):
                # Nested dict: process its key-values under 'pos'
                for subkey, subvalue in value.items():
                    predicate = URIRef(f"{fb_namespace}{subkey}")
                    self._add_literal(pos_uri, predicate, subvalue, datatype=XSD.decimal)
            else:
                # All other values belong directly under 'position'
                predicate = URIRef(f"{fb_namespace}{key}")
                self._add_literal(position_uri, predicate, value)
    
    def _process_level(self, parent_uri, level_data, fx_namespace, fb_namespace):
        """
        Create RDF triples for a structured flight level object.
        Handles cases where flightLevel includes both value and uom.
        
        Expected level_data format:
        {
            'flightLevel': {
                'uom': 'FL',
                'flightLevel': '370'  # or int
            }
        }
        """
        level_uri = URIRef(f"{parent_uri}_level")

        # Link parent to level
        self.graph.add((parent_uri, fx_namespace.level, level_uri))
        self.graph.add((level_uri, RDF.type, fx_namespace.Level))

        for key, value in level_data.items():
            if key == "flightLevel" and isinstance(value, dict):
                # Handle both value and uom
                for subkey, subval in value.items():
                    pred = URIRef(f"{fb_namespace}{subkey}")
                    if subkey == "flightLevel":
                        self._add_literal(level_uri, pred, subval, datatype=XSD.integer)
                    else:
                        self._add_literal(level_uri, pred, subval)
            else:
                # Handle any unexpected structure
                pred = URIRef(f"{fb_namespace}{key}")
                self._add_literal(level_uri, pred, value)


    def _process_time(self, parent_uri, time_data, fx_namespace, fb_namespace):
        """
        Create RDF triples for a structured flight time object.
        - Output structure: parent -> time
        """
        time_uri = URIRef(f"{parent_uri}_time")

        # Structure: parent -> time
        self.graph.add((parent_uri, fx_namespace.time, time_uri))
        self.graph.add((time_uri, RDF.type, fx_namespace.Time))

        for key, value in time_data.items():
            predicate = URIRef(f"{fb_namespace}{key}")
            self._add_literal(time_uri, predicate, value, datatype=XSD.dateTime) # dateTime

    def _process_hmi_position(self, flight_uri, uri, key, value, namespace):
        """
        Create RDF triples for a structured HMI position object.
        - Keys with nested dicts are processed under 'pos'
        - All other top-level keys are added under 'position'
        """

        self.graph.add((uri, RDF.type, HMI.Label))
        parent_predicate = URIRef(HMI[key])
        self.graph.add((flight_uri, parent_predicate, uri))

        for subkey, subvalue in value.items():
            if subvalue:
                predicate = URIRef(f"{namespace}{subkey}")
                self._add_literal(uri, predicate, subvalue, datatype=XSD.integer)
    
    def _process_popup(self, flight_uri, uri, key, value, namespace):
        """
        Create RDF triples for a structured HMI popup object.
        - Keys with nested dicts are processed under 'popup'
        - All other top-level keys are added under 'popup'
        """

        self.graph.add((uri, RDF.type, HMI.Popup))
        parent_predicate = URIRef(HMI[key])
        self.graph.add((flight_uri, parent_predicate, uri))

        for subkey, subvalue in value.items():
            if subvalue:
                predicate = URIRef(f"{namespace}{subkey}")
                self._add_literal(uri, predicate, subvalue)

    # TODO: merge popup and hmi position processing methods into one

    def _process_mouse_position(self, hmi_uri, mouse_uri, hmi_entity_data):
        self.graph.add((mouse_uri, RDF.type, HMI.MousePosition))
        self.graph.add((hmi_uri, HMI.mousePosition, mouse_uri))
        for key, value in hmi_entity_data.items():
            self._add_literal(mouse_uri, URIRef(HMI[key]), value, datatype=XSD.integer)

    def _process_alert(self, hmi_uri, alert_uri, hmi_entity_data):
        self.graph.add((alert_uri, RDF.type, HMI.Alert))
        self.graph.add((hmi_uri, HMI.alert, alert_uri))

        for key, value in hmi_entity_data.items():
            if key == "relatedFlights" and isinstance(value, dict):
                for _, flight_id in value.items():
                    flight_uri = URIRef(FIXM[flight_id])
                    self.graph.add((alert_uri, HMI.relatedFlight, flight_uri))
                    self.graph.add((flight_uri, RDF.type, FIXM.Flight))
            else:
                self._add_literal(alert_uri, URIRef(HMI[key]), value)

    def _process_velocity_vector(self, point4d_uri, velocity_data):
        """Process velocity data and add to the graph."""
        velocity_uri = URIRef(f"{point4d_uri}_predictedAirspeed")
        self.graph.add((velocity_uri, RDF.type, FIXM.predictedAirspeed))
        self.graph.add((point4d_uri, FIXM.predictedAirspeed, velocity_uri))
        self._add_literal(velocity_uri, FIXM.Vx, velocity_data.get("VxMs"), datatype=XSD.decimal)
        self._add_literal(velocity_uri, FIXM.Vy, velocity_data.get("VyMs"), datatype=XSD.decimal)

    def _process_clearance(self, point4d_uri, point4d_data):
        """Process clearance info and add to the graph."""
        print("point4d_Data: ", point4d_data)
        for key, value in point4d_data.items():
            print("Key in clearance: ", key)
            if key in self.clearance_keys: # TODO: add other clearance types
                clearance_uri = URIRef(f"{point4d_uri}_{key}")
                self.graph.add((clearance_uri, RDF.type, URIRef(FIXM[key])))
                self.graph.add((point4d_uri, URIRef(FIXM[key]), clearance_uri))
                self._add_literal(clearance_uri, URIRef(FIXM[key]), value)

    def _process_point4d(self, element_uri, point4d_data):
        """Process point4D data and add to the graph."""
        point4d_uri = URIRef(f"{element_uri}_point4D")
        self.graph.add((element_uri, FIXM.point4D, point4d_uri))
        self.graph.add((point4d_uri, RDF.type, FIXM.Point4D))

        print("Point4d data: ", point4d_data)

        for key, value in point4d_data.items():
            if key == "position":
                self._process_position(point4d_uri, value, FIXM, FB)

            elif key == "level":
                self._process_level(point4d_uri, value, FIXM, FB)

            elif key == "time":
                self._process_time(point4d_uri, value, FIXM, FB)

            elif key == "velocity":
                self._process_velocity_vector(point4d_uri, value)

            elif key == "flightLevelM":
                self._add_literal(point4d_uri, FIXM.flightLevel, value, datatype=XSD.decimal)

            elif key == "calculatedAltitudeM":
                self._add_literal(point4d_uri, FIXM.calculatedAltitude, value, datatype=XSD.decimal)
            
            else:
                # For all other keys, process as clearance or custom property
                self._process_clearance(point4d_uri, {key: value})
    
    def _process_branch(self, route_trajectory_group_uri, branch_key, branch_data):
        branch_uri = URIRef(f"{route_trajectory_group_uri}_{branch_key}")
        self.graph.add((route_trajectory_group_uri, URIRef(FIXM[branch_key]), branch_uri))
        self.graph.add((branch_uri, RDF.type, URIRef(FIXM[branch_key.capitalize()])))

        if branch_key in ["agreed", "desired"]:
            branch_data = branch_data[0] # Take only the first point as the next point

        element = branch_data.get("element", {})

        element_uri = URIRef(f"{branch_uri}_element")
        self.graph.add((branch_uri, FIXM.element, element_uri))
        self.graph.add((element_uri, RDF.type, FIXM.Element))

        designator = (
            element
            .get("elementStartPoint", {})
            .get("designatedPoint", {})
            .get("designator", None)
        )
        if designator:
            self.graph.add((element_uri, FIXM.designator, Literal(designator)))

        point4d = element.get("point4D", {})
        self._process_point4d(element_uri, point4d)

    def _process_predicted_trajectory(self, route_trajectory_group_uri, trajectory_data):
        """Process trajectory data and add to the graph."""
        print("trajectory_data: ", trajectory_data)
        predicted_uri = URIRef(f"{route_trajectory_group_uri}_predicted")
        self.graph.add((route_trajectory_group_uri, FIXM.predicted, predicted_uri))
        self.graph.add((predicted_uri, RDF.type, FIXM.Predicted))

        # TODO: prefix shouldn't be fx here - update this
        self._add_literal(predicted_uri, FIXM.points, trajectory_data)

    def _process_label(self, flight_uri, label_data):
        """Process label data and add to the graph."""
        label_uri = URIRef(f"{flight_uri}_label")
        self.graph.add((label_uri, RDF.type, HMI.Label))
        self.graph.add((flight_uri, HMI.label, label_uri))

        if "position" in label_data:
            self._process_position(label_uri, label_data["position"], HMI)
    
    def _process_route_trajectory_group(self, flight_uri, value):
        """Process route trajectory group data and add to the graph."""
        route_trajectory_group_uri = URIRef(f"{flight_uri}_RTG")
        self.graph.add((flight_uri, FIXM.routeTrajectoryGroup, route_trajectory_group_uri))
        self.graph.add((route_trajectory_group_uri, RDF.type, FIXM.RouteTrajectoryGroup))

        for branch_key, branch_data in value.items():
            if branch_key == "predicted":
                self._process_predicted_trajectory(route_trajectory_group_uri, branch_data)
            else:
                self._process_branch(route_trajectory_group_uri, branch_key, branch_data)

    def _process_non_route_data(self, uri, key, value):
        if isinstance(value, dict):
            # If value is a dictionary, create a new URI and process it recursively
            sub_uri = URIRef(f"{uri}_{key}")
            self.graph.add((uri, URIRef(FIXM[key]), sub_uri))
            for sub_key, sub_value in value.items():
                self._process_non_route_data(sub_uri, sub_key, sub_value)
        else:
            # If value is a literal, add it directly
            self._add_literal(uri, URIRef(FIXM[key]), value)

    def _process_flight(self, flight_uri, flight_data):
        """Process flight data and add to the graph."""

        for key, value in flight_data.items():
            # Main flight-related data processing
            if key == "routeTrajectoryGroup":
                self._process_route_trajectory_group(flight_uri, value)
            
            # Handle other flight-related data
            else:
                self._process_non_route_data(flight_uri, key, value)

    def _process_hmi(self, timestamp_uri, hmi_data):
        """Process HMI data and add to the graph."""

        hmi_uri = URIRef(f"{timestamp_uri}_hmi")
        self.graph.add((hmi_uri, RDF.type, HMI.Interface))
        self.graph.add((timestamp_uri, HMI.interface, hmi_uri))

        for hmi_entity_key, hmi_entity_data in hmi_data.items():
            if hmi_entity_key.lower().startswith("flight_"):
                flight_id = hmi_entity_key.split("_", 1)[1]
                flight_uri = URIRef(f"{hmi_uri}_flight_{flight_id}")
                self.graph.add((hmi_uri, HMI.relatedFlight, flight_uri))
                self.graph.add((flight_uri, RDF.type, HMI.FlightRepresentation))

                for key, value in hmi_entity_data.items():
                    if key == "trackLabelPosition":
                        label_uri = URIRef(f"{flight_uri}_labelPos")
                        self._process_hmi_position(flight_uri, label_uri, key, value, HMI)

                    elif key == "trackScreenPosition":
                        screen_uri = URIRef(f"{flight_uri}_screenPos")
                        self._process_hmi_position(flight_uri, screen_uri, key, value, HMI)

                    elif key == "popup":
                        popup_uri = URIRef(f"{flight_uri}_popup")
                        self._process_popup(flight_uri, popup_uri, key, value, HMI)

                    else:
                        # TODO: does this case exist?
                        self._add_literal(flight_uri, URIRef(HMI[key]), value)

            elif hmi_entity_key == "mousePosition":
                mouse_uri = URIRef(f"{hmi_uri}_mousePosition")
                self._process_mouse_position(hmi_uri, mouse_uri, hmi_entity_data)

            elif hmi_entity_key == "alert":
                alert_uri = URIRef(f"{hmi_uri}_alert")
                self._process_alert(hmi_uri, alert_uri, hmi_entity_data)

    def _repository_to_rdf(self, timestamp_data):
        """
        Converts the data repository into an RDF graph.
        This method takes the most recent timestamp from a data repository and converts data into RDF triples.
        """
        timestamp, data = next(iter(timestamp_data.items()))
        timestamp = self._safe_local_name(timestamp)
        timestamp_uri = BASE[f"{timestamp}"]
        self.graph.add((timestamp_uri, RDF.type, BASE.Timestamp))

        for key, key_data in data.items():
            if key.startswith("flight_"):
                flight_uri = URIRef(FIXM[f"{key}_{timestamp}"])
                self.graph.add((timestamp_uri, FIXM.flight, flight_uri))
                self.graph.add((flight_uri, RDF.type, FIXM.Flight))
                if "Flight" in key_data:
                    self._process_flight(flight_uri, key_data["Flight"])

            elif key == "HMI":
                # Process HMI data
                hmi_data = key_data
                self._process_hmi(timestamp_uri, hmi_data)
    

    #########################################################################################################
    ###################################     DATA REPOSITORY METHODS     #####################################
    #########################################################################################################
    
    
    def _safe_local_name(self, s):
        return s.replace("-", "_").replace(":", "_").replace(".", "_")

    def _match_uuids_and_callsigns(self, uuid: str, xml_data: str):
        """
        Parses XML data and maps flight UUID to its corresponding aircraftIdentification (callsign).

        Args:
            uuid (str): UUID of the flight plan to match.
            xml_data (str): String containing the full Flight Plan XML.
        """

        # Convert FIXM Namespace to string URI (ElementTree requires raw string URIs, not rdflib Namespace objects)
        fixm_uri = str(self.namespaces["fx"])
        fixm_uri = fixm_uri.rstrip("/")  # Ensure no trailing slash

        # Parse XML
        root = ET.fromstring(xml_data)

        # Use fully qualified tag name to locate the element
        element = root.find(f".//{{{fixm_uri}}}aircraftIdentification")
        
        if element is not None and element.text:
            callsign = element.text.strip()
            self.uuid_callsign_map[uuid] = callsign
        else:
            print("Warning: aircraftIdentification element not found.")
            callsign = None

        self.aircraft_id = callsign  # Store the aircraft ID for later use

    def _update_asd_event_repo(self, json_record):
        data_record = json_record.get("asd_event", {})
        timestamp = data_record.get("timestamp")  # TODO: which timestamp should we use here? asd_event or the one from the main record?

        # 1. Handle clearance 
        clearance = data_record.get("clearance")
        if clearance:
            track_number = clearance.get("trackNumber")
            if track_number:
                callsign = self.track_number_callsign_map.get(track_number)
                flight_key = f"flight_{callsign}"

                clearance_type = clearance.get("clearanceType")
                clearance_value = clearance.get("clearance")

                if callsign and clearance_type and clearance_value:
                    # Navigate safely to "agreed" list
                    route_group = self.data_repository \
                        .setdefault(timestamp, {}) \
                        .setdefault(flight_key, {}) \
                        .setdefault("Flight", {}) \
                        .setdefault("routeTrajectoryGroup", {})

                    agreed_list = route_group.get("agreed")
                    if not isinstance(agreed_list, list):
                        print(f"Warning: 'agreed' field is not a list. Skipping.")
                        return

                    if not agreed_list:
                        # If the list is empty, add a new element structure
                        agreed_list.append({"element": {"point4D": {}}})

                    # Use first element in the list for now
                    first_element = agreed_list[0].get("element", {})
                    point4d = first_element.setdefault("point4D", {})

                    # Update clearance
                    # TODO: Make sure this updates exactly a directory/subdirectory it refers to
                    point4d[clearance_type] = clearance_value

                    # Write back the updated element in case it was not present before
                    agreed_list[0]["element"] = first_element

        # 2. Determine if there's a flight-related HMI event
        track_number = (
            data_record.get("trackLabelPosition", {}).get("trackNumber")
            or data_record.get("trackScreenPosition", {}).get("trackNumber")
            or data_record.get("popup", {}).get("trackNumber")
        )

        if track_number:
            callsign = self.track_number_callsign_map.get(track_number)
            if callsign:
                flight_key = f"flight_{callsign}"
                hmi_branch = self.data_repository \
                    .setdefault(timestamp, {}) \
                    .setdefault("HMI", {}) \
                    .setdefault(flight_key, {})

                hmi_branch["aircraftIdentification"] = callsign

                # 2a. trackLabelPosition
                if "trackLabelPosition" in data_record:
                    tlp = data_record["trackLabelPosition"]
                    hmi_branch["trackLabelPosition"] = {
                        "x": tlp.get("x"),
                        "y": tlp.get("y"),
                        "width": tlp.get("width"),
                        "height": tlp.get("height"),
                        "visible": tlp.get("visible"),
                        "selected": tlp.get("selected"),
                    }

                # 2b. trackScreenPosition
                if "trackScreenPosition" in data_record:
                    tsp = data_record["trackScreenPosition"]
                    hmi_branch["trackScreenPosition"] = {
                        "x": tsp.get("x"),
                        "y": tsp.get("y"),
                    }

                # 2c. popup
                if "popup" in data_record:
                    popup = data_record["popup"]
                    hmi_branch["popup"] = {
                        "name": popup.get("name"),
                        "opened": popup.get("opened"),
                    }

        # 3. Handle mousePosition (not tied to specific flight)
        if "mousePosition" in data_record:
            mouse_position = data_record["mousePosition"]
            self.data_repository \
                .setdefault(timestamp, {}) \
                .setdefault("HMI", {})["mousePosition"] = mouse_position

    def _update_alert_repo(self, json_record):
        data_record = json_record.get("alert", {})
        timestamp = json_record.get("timestamp")

        # Get track numbers and their corresponding callsigns
        track_numbers = [data_record[key] for key in data_record if "trackNumber" in key]
        related_callsigns = [self.track_number_callsign_map[track_num] for track_num in track_numbers]

        # Ensure alert structure exists in the repository
        alert_group = self.data_repository \
            .setdefault(timestamp, {}) \
            .setdefault("HMI", {}) \
            .setdefault("alert", {})

        # Update relatedFlights (indexed callsign references)
        alert_group.setdefault("relatedFlights", {})
        for idx, callsign in enumerate(related_callsigns):
            alert_group["relatedFlights"][idx] = f"Flight_{callsign}"

        # Add other alert fields (excluding trackNumbers)
        for key, value in data_record.items():
            if "trackNumber" not in key:
                alert_group[key] = value

    def _update_track_repo(self, json_record):
        timestamp = json_record.get("timestamp")
        data_record = json_record.get("track", {})
        callsign = data_record.get("callSign")
        flight_key = f"flight_{callsign}"

        # Build point4D using nested structure from track record
        point4d = {}

        # Position
        position = data_record.get("position")
        if position:
            lat = position.get("latDeg")
            lon = position.get("lonDeg")
            if lat is not None and lon is not None:
                point4d["position"] = {
                    "srsName": "urn:ogc:def:crs:EPSG::4326",
                    "pos": {
                        "lat": lat,
                        "lon": lon
                    }
                }

        # Velocity
        velocity = data_record.get("velocity")
        if velocity:
            vx = velocity.get("VxMs")
            vy = velocity.get("VyMs")
            if vx is not None or vy is not None:
                point4d["velocity"] = {}
                if vx is not None:
                    point4d["velocity"]["VxMs"] = vx
                if vy is not None:
                    point4d["velocity"]["VyMs"] = vy

        # Level block (with FL values)
        flight_level_m = data_record.get("flightLevelM")
        calculated_alt_m = data_record.get("calculatedAltitudeM")

        level = {}
        if flight_level_m is not None:
            level["flightLevel"] = str(int(round(flight_level_m * self.to_fl_const)))
        if calculated_alt_m is not None:
            level["calculatedFlightLevel"] = str(int(round(calculated_alt_m * self.to_fl_const)))

        if level:
            point4d["level"] = level

        # Final repo update
        route_group = self.data_repository \
            .setdefault(timestamp, {}) \
            .setdefault(flight_key, {}) \
            .setdefault("Flight", {}) \
            .setdefault("routeTrajectoryGroup", {})

        route_group["current"] = {
            "element": {
                "point4D": point4d
            }
        }

    def _update_predicted_trajectory_repo(self, json_record):
        timestamp = json_record.get("timestamp")
        data_record = json_record.get("trajectory", {})
        flight_plan_uuid = data_record.get("flightPlanUuid")
        callsign = self.uuid_callsign_map.get(flight_plan_uuid)
        
        if not callsign:
            print(f"Warning: No callsign found for UUID {flight_plan_uuid}. Skipping predicted trajectory update.")
            return  # or raise a warning/log error if UUID not found

        flight_key = f"flight_{callsign}"
        predicted_points = data_record.get("points", [])

        # Ensure the structure exists and only update the "predicted" section
        route_group = self.data_repository \
            .setdefault(timestamp, {}) \
            .setdefault(flight_key, {}) \
            .setdefault("Flight", {}) \
            .setdefault("routeTrajectoryGroup", {})

        route_group["predicted"] = predicted_points
        print("Route group after predicted update: ", route_group)
        
    def _update_data_repository_events(self, json_record):
        """
        Processes an aviation-related event record data and organizes it into a structured dictionary by timestamp and flight.

        The function performs the following steps:
        1. Builds a mapping from track numbers to callsigns by scanning all track events.
        2. Extract and organize information such as flight tracks, clearances, trajectories, HMI (Human-Machine Interface) events, and alerts.
        3. For each timestamp, it creates a nested dictionary structure containing flight-specific and HMI-related data.
        4. For records referencing a flight by UUID, it uses the provided dictionary to resolve the callsign.
        5. Returns both the structured dictionary and the track-to-callsign mapping.

        Args:
            json_record (dict/json): Json data representing an aviation event or record.
        Returns:
            tuple: (timestamp_repository, track_callsign_dict)
                - timestamp_repository (dict): Nested dictionary organized by timestamp and flight/HMI data.
        """
        # Extract timestamp from appropriate location in record
        if "asd_event" in json_record:
            timestamp = json_record.get("asd_event").get("timestamp")
        else:
            timestamp = json_record.get("timestamp")

        # Copy previous repo state to the current timestamp
        self._create_new_repo_state(timestamp)
        # Remove copied predicted trajectories from the previous timestamp to save memmory
        self._remove_predicted_trajectories_from_last_state_repo(timestamp)

        # Process track events to fill the track-callsign mapping dictionary
        track = json_record.get("track", {})
        track_number = track.get("trackNumber")
        call_sign = track.get("callSign")

        if track_number is not None and call_sign is not None:
            self.track_number_callsign_map.update(
                {track_number: call_sign}
            )

        if "asd_event" in json_record:
            # Process asd_event records
            self._update_asd_event_repo(json_record)
                
        if "alert" in json_record:
            # Process alert records
            self._update_alert_repo(json_record)

        if "track" in json_record:
            # Process track/current position records
            self._update_track_repo(json_record)
            
        elif "trajectory" in json_record:
            # Process predicted trajectory records
            self._update_predicted_trajectory_repo(json_record)

        self.last_timestamp = timestamp  # Store the last processed timestamp

        return self.data_repository, self.track_number_callsign_map
    
    def _update_data_repository_fixm(self, fixm_record, timestamp, callsign):
        """
        Parses a FIXM-compliant XML flight record and converts it into a clean, nested dictionary
        with the following transformations:
            - XML attributes are stored as regular keys.
            - Text values are stored directly.
            - Nested child elements are recursively processed.
            - Repeated tags are converted to lists.
        
        The resulting dictionary is stored under a key named `flight_<callsign>`.
        """

        # Copy previous repo state to the current timestamp
        self._create_new_repo_state(timestamp)
        # Remove copied predicted trajectories from the previous timestamp to save memory
        self._remove_predicted_trajectories_from_last_state_repo(timestamp)

        try:
            root = ET.fromstring(fixm_record)
        except ET.ParseError as e:
            print(f"XML Parsing Error in _update_data_repository_fixm: {e}")
            return {}

        def element_to_dict(element, parent_tag=None):
            tag = element.tag.split("}")[-1]
            children = list(element)

            if tag in {"agreed", "desired"}:
                node = []
                for child in children:
                    child_tag = child.tag.split("}")[-1]
                    child_dict = element_to_dict(child, parent_tag=tag)
                    if child_tag == "element":
                        node.append({child_tag: child_dict})
                return node

            # Normal case
            node = {}

            # Store attributes first
            for attr_name, attr_value in element.attrib.items():
                node[attr_name] = attr_value

            if children:
                for child in children:
                    child_tag = child.tag.split("}")[-1]
                    child_dict = element_to_dict(child, parent_tag=tag)

                    if child_tag not in node:
                        node[child_tag] = child_dict
                    else:
                        if not isinstance(node[child_tag], list):
                            node[child_tag] = [node[child_tag]]
                        node[child_tag].append(child_dict)
            else:
                text = element.text.strip() if element.text else ""

                if tag == "pos" and text:
                    try:
                        lat_str, lon_str = text.split()
                        return {"lat": float(lat_str), "lon": float(lon_str)}
                    except ValueError:
                        return text

                # ⬇️ Merge text with attributes if both exist
                if text:
                    if node:  # There were attributes
                        node[tag] = text
                        return node
                    else:
                        return text
                elif node:
                    return node
                else:
                    return ""

            return node

        # Top-level dictionary with root tag
        root_tag = root.tag.split("}")[-1]
        structured_dict = {root_tag: element_to_dict(root)}

        # Store under a timestamp-based key or unique flight key
        flight_key = f"flight_{callsign}"

        if flight_key in self.data_repository[timestamp].keys():
            self.data_repository[timestamp][flight_key].update(structured_dict)
        else:
            self.data_repository[timestamp].update({
                flight_key: structured_dict
            })

        self.last_timestamp = timestamp


    def _create_new_repo_state(self, new_timestamp):
        """Create a new repository state for the next timestamp."""
        previous_state = self.data_repository.get(self.last_timestamp, {})
        self.data_repository[new_timestamp] = copy.deepcopy(previous_state)

        # print(f"Created new repository state for timestamp: {new_timestamp}, used previous state: {self.last_timestamp}")

    def _remove_predicted_trajectories_from_last_state_repo(self, new_timestamp):
        """
        Avoid copying predicted trajectories from the previous timestamp to the new one.
        This is to avoid copying large strings - we can look them up in history timestamps if needed.
        This method modifies the last state in the data repository by removing the "predicted" key.
        """

        # TODO: optimise this not to check all flights (for loop below) for every new timestamp
        # - use some flag to check if predicted trajectories were added and processed

        print(f"Removing predicted trajectories from last state for timestamp: {new_timestamp}")

        if new_timestamp in self.data_repository:
            last_state = self.data_repository[new_timestamp]
            print(f"Last state for timestamp {new_timestamp}: {last_state}")
            for flight_key, flight_data in last_state.items():
                if "Flight" in flight_data and "routeTrajectoryGroup" in flight_data["Flight"]:
                    route_group = flight_data["Flight"]["routeTrajectoryGroup"]
                    # Remove predicted trajectory from the new state if exists
                    route_group.pop("predicted", None)
        
    #########################################################################################################
    ##########################################     MAIN METHODS     #########################################
    #########################################################################################################

    def convert_event_data(self, json_record):
        """Convert non-FIXM data to RDF graph."""

        # This handles the data repository creation
        self._update_data_repository_events(json_record) # TODO: can we convert repository to turtle without additional rdf graph creation?
        # print(json.dumps(self.data_repository, indent=4))
                
    def convert_fixm_data(self, xml_record):
        """Convert FIXM/XML data to RDF graph."""
        timestamp = xml_record["timestamp"]
        flight_plan = xml_record["flight_plan"]
        fixm_data = flight_plan["fixm"]
        uuid = flight_plan["uuid"]
        self._match_uuids_and_callsigns(uuid, fixm_data)
        matched_callsign = self.uuid_callsign_map.get(uuid)
        
        # This handles the data repository creation
        self._update_data_repository_fixm(fixm_data, timestamp, matched_callsign) # TODO: can we convert repository to turtle without additional rdf graph creation?
        # print(json.dumps(self.data_repository, indent=4))

    def convert_record_data(self, data_record, data_record_type):
        """Convert a record (either JSON or XML) to RDF graph."""

        if data_record_type == "json":
            # 1. case - FIXM data
            if "flight_plan" in data_record.keys():
                self.convert_fixm_data(data_record)
                # TODO: check what happens if there are two same timestamps for xml and json, will they overwrite each other?
            
            # 2. case - Event data
            else:
                self.convert_event_data(data_record)
                # TODO: check what happens if there are two same timestamps for xml and json, will they overwrite each other?
        
            timestamp_data = {self.last_timestamp: self.data_repository[self.last_timestamp]}
            # Convert the repository to RDF triples
            self._repository_to_rdf(timestamp_data)
        else:
            raise ValueError("Unsupported data record type. Use 'json'.")

    def serialize(self, format="turtle"):
        """Serialize the RDF graph to a string/turtle format."""
        return self.graph.serialize(format=format)

    # For testing purposes, you can save the RDF graph to a file
    def save(self, file_path, format="turtle"):
        """Save the RDF graph to a turtle file."""
        with open(file_path, "w", encoding="utf-8") as out_file:
            out_file.write(self.serialize(format=format))