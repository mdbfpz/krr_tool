from rdflib import Graph, URIRef, Literal, Namespace
from rdflib.namespace import RDF, RDFS, XSD 
import xml.etree.ElementTree as ET
import json
import copy

# Define namespaces
FIXM = Namespace("http://www.fixm.aero/flight/4.3")
FB = Namespace("http://www.fixm.aero/base/4.3")
XS = Namespace("http://www.w3.org/2001/XMLSchema")
HMI = Namespace("http://www.fixm.aero/hmi/1.0/") # TODO: this should not be a part of FIXM, but a separate namespace, maybe 'plain' extension

class RDFConverter:
    def __init__(self):
        """Initialize the RDF graph and bind namespaces."""
        self.graph = Graph()
        self.graph.bind("fx", FIXM)
        self.graph.bind("fb", FB)
        self.graph.bind("hmi", HMI)
        self.graph.bind("xs", XS)
        self.graph.bind("xsd", XSD)
        self.graph.bind("rdf", RDF)
        self.graph.bind("rdfs", RDFS)

        # TODO: remove this if possible and use variables declared above    
        self.fixm_namespaces = {
            "fx": FIXM,
            "fb": FB,
        }

        self.data_repository = {}
        self.track_number_callsign_map = {} # key = track number; value = callsign
        self.uuid_callsign_map = {}  # key = flight plan UUID; value = callsign
        self.last_timestamp = None  # Store the last processed timestamp

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

    def _add_uri(self, subject, predicate, obj):
        """Add a URI triple to the graph."""
        self.graph.add((subject, predicate, obj))

    def _add_triples(self, subject, predicates):
        for predicate_uri, objects in predicates.items():
            predicate = URIRef(predicate_uri)
            for obj in objects:
                value = obj["value"]
                if value == "literal":
                    self._add_literal(subject, predicate, value)
                elif value == "uri":
                    obj_uri = URIRef(value)
                    self._add_uri(subject, predicate, obj_uri)

    def _process_position(self, parent_uri, position_data, namespace):
        """Process position data and add to the graph."""
        position_uri = URIRef(f"{parent_uri}/position")
        self.graph.add((position_uri, RDF.type, namespace.Position))
        self.graph.add((parent_uri, namespace.position, position_uri))
        for key, value in position_data.items():
            self._add_literal(position_uri, URIRef(f"{namespace}{key}"), value)

    def _process_velocity(self, point4d_uri, velocity_data):
        """Process velocity data and add to the graph."""
        velocity_uri = URIRef(f"{point4d_uri}/predictedAirspeed")
        self.graph.add((velocity_uri, RDF.type, FIXM.predictedAirspeed))
        self.graph.add((point4d_uri, FIXM.predictedAirspeed, velocity_uri))
        self._add_literal(velocity_uri, FIXM.Vx, velocity_data.get("vx_ms"), XSD.decimal)
        self._add_literal(velocity_uri, FIXM.Vy, velocity_data.get("vy_ms"), XSD.decimal)

    def _process_point4d(self, element_uri, point4d_data):
        """Process point4D data and add to the graph."""
        point4d_uri = URIRef(f"{element_uri}/point4D")
        self.graph.add((element_uri, FIXM.point4D, point4d_uri))
        self.graph.add((point4d_uri, RDF.type, FIXM.Point4D))

        if "position" in point4d_data:
            self._process_position(point4d_uri, point4d_data["position"], FIXM)

        if "velocity" in point4d_data:
            self._process_velocity(point4d_uri, point4d_data["velocity"])

        if "flight_level_m" in point4d_data:
            self._add_literal(point4d_uri, FIXM.flightLevel, point4d_data["flight_level_m"], XSD.decimal)

        if "calculated_altitude_m" in point4d_data:
            self._add_literal(point4d_uri, FIXM.calculatedAltitude, point4d_data["calculated_altitude_m"], XSD.decimal)

    def _process_current(self, flight_uri, current_data):
        """Process current flight data and add to the graph."""
        current_uri = URIRef(f"{flight_uri}/current")
        self.graph.add((flight_uri, FIXM.current, current_uri))
        self.graph.add((current_uri, RDF.type, FIXM.Current))

        if "element" in current_data:
            element_uri = URIRef(f"{current_uri}/element")
            self.graph.add((current_uri, FIXM.element, element_uri))
            self.graph.add((element_uri, RDF.type, FIXM.Element))

            if "point4D" in current_data["element"]:
                self._process_point4d(element_uri, current_data["element"]["point4D"])

    def _process_trajectory(self, flight_uri, trajectory_data):
        """Process trajectory data and add to the graph."""
        trajectory_uri = URIRef(f"{flight_uri}/trajectory")
        self.graph.add((flight_uri, FIXM.predicted, trajectory_uri))
        self.graph.add((trajectory_uri, RDF.type, FIXM.Trajectory))

        if "points" in trajectory_data:
            lat_lon_fl_str = " ".join(
                f"{point['lat_deg']} {point['lon_deg']} {point['flight_level']}" for point in trajectory_data["points"]
            )
            self._add_literal(trajectory_uri, FIXM.trajectoryPoints, lat_lon_fl_str)

    def _process_label(self, flight_uri, label_data):
        """Process label data and add to the graph."""
        label_uri = URIRef(f"{flight_uri}/label")
        self.graph.add((label_uri, RDF.type, HMI.Label))
        self.graph.add((flight_uri, HMI.label, label_uri))

        if "position" in label_data:
            self._process_position(label_uri, label_data["position"], HMI)

    def _process_flight(self, timestamp_uri, entity_key, entity_data):
        """Process flight data and add to the graph."""
        flight_id = entity_key.split("_")[1]
        flight_uri = URIRef(FIXM[f"flight_{flight_id}"])
        self.graph.add((timestamp_uri, FIXM.flight, flight_uri))
        self.graph.add((flight_uri, RDF.type, FIXM.Flight))

        if "flightIdentification" in entity_data:
            self._add_literal(flight_uri, FIXM.aircraftIdentification, entity_data["flightIdentification"])

        if "current" in entity_data:
            self._process_current(flight_uri, entity_data["current"])

        if "trajectory" in entity_data:
            self._process_trajectory(flight_uri, entity_data["trajectory"])

        if "label" in entity_data:
            self._process_label(flight_uri, entity_data["label"])

    def _process_hmi(self, timestamp_uri, hmi_data):
        """Process HMI data and add to the graph."""
        hmi_uri = URIRef(f"{timestamp_uri}/hmi")
        self.graph.add((hmi_uri, RDF.type, HMI.Interface))
        self.graph.add((timestamp_uri, HMI.interface, hmi_uri))

        for hmi_entity_key, hmi_entity_data in hmi_data.items():
            if hmi_entity_key.startswith("Flight_"):
                flight_id = hmi_entity_key.split("_")[1]
                flight_uri = URIRef(f"hmi/flight_{flight_id}")
                self.graph.add((hmi_uri, HMI.relatedFlight, flight_uri))

                if "track_label_position" in hmi_entity_data:
                    label_uri = URIRef(f"{flight_uri}/label")
                    self.graph.add((label_uri, RDF.type, HMI.Label))
                    if "position" in hmi_entity_data["track_label_position"]:
                        self._process_position(label_uri, hmi_entity_data["track_label_position"]["position"], HMI)

                if "track_screen_position" in hmi_entity_data:
                    label_uri = URIRef(f"{flight_uri}/screen")
                    self.graph.add((label_uri, RDF.type, HMI.Label))
                    if "position" in hmi_entity_data["track_screen_position"]:
                        self._process_position(label_uri, hmi_entity_data["track_screen_position"]["position"], HMI)

            elif hmi_entity_key == "mouse_position":
                mouse_uri = URIRef(f"{hmi_uri}/mousePosition")
                self.graph.add((mouse_uri, RDF.type, HMI.MousePosition))
                self.graph.add((hmi_uri, HMI.mousePosition, mouse_uri))
                for key, value in hmi_entity_data.items():
                    self._add_literal(mouse_uri, URIRef(f"{HMI}{key}"), value)

            elif hmi_entity_key == "alert":
                alert_uri = URIRef(f"{hmi_uri}/alert")
                self.graph.add((alert_uri, RDF.type, HMI.alert))
                self.graph.add((hmi_uri, HMI.alert, alert_uri))

                for key, value in hmi_entity_data.items():
                    if key == "relatedFlights":
                        related_flights_uri = URIRef(f"{alert_uri}/relatedFlights")
                        for flight_key, flight_value in value.items():
                            flight_uri = URIRef(FIXM[f"{flight_value}"])
                            self.graph.add((related_flights_uri, RDF.type, flight_uri))
                            self.graph.add((flight_uri, RDF.type, FIXM.Flight))
                    else:
                        self._add_literal(alert_uri, URIRef(f"{HMI}{key}"), value)

    def _match_uuids_and_callsigns(self, uuid: str, xml_data: str) -> dict:
        """
        Parses XML data and maps flight UUID to its corresponding aircraftIdentification (callsign).

        Args:
            uuid (str): UUID of the flight plan to match.
            xml_data (str): String containing the full Flight PLan XML.

        Returns:
            dict: Mapping from UUID to aircraftIdentification value.
        """

        root = ET.fromstring(xml_data)
        callsign = root.find(
            ".//fx:aircraftIdentification", self.fixm_namespaces
        ).text
        print("Parsed callsign: ", callsign)
        self.uuid_callsign_map[uuid] = callsign

        self.aircraft_id = callsign  # Store the aircraft ID for later use
    
    def _update_asd_event_repo(self, json_record):
        data_record = json_record.get("asd_event", {})
        timestamp = data_record.get("timestamp")  # TODO: which timestamp should we use here? asd_event or the one from the main record?

        # === 1. Handle clearance ===
        # TODO: test if this works
        clearance = data_record.get("clearance")
        if clearance:
            track_number = clearance.get("trackNumber")
            if track_number:
                callsign = self.track_number_callsign_map.get(track_number)
                flight_key = f"flight_{callsign}"

                # Set clearance value under agreed → element → point4D → [type]
                clearance_type = clearance.get("clearanceType")
                clearance_value = clearance.get("clearance")
                if callsign and clearance_type and clearance_value:
                    flight_branch = self.data_repository \
                        .setdefault(timestamp, {}) \
                        .setdefault(flight_key, {}) \
                        .setdefault("agreed", {}) \
                        .setdefault("element", {}) \
                        .setdefault("point4D", {})

                    flight_branch[clearance_type] = clearance_value

        # === 2. Determine if there's a flight-related HMI event ===
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

                # === 2a. trackLabelPosition ===
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

                # === 2b. trackScreenPosition ===
                if "trackScreenPosition" in data_record:
                    tsp = data_record["trackScreenPosition"]
                    hmi_branch["trackScreenPosition"] = {
                        "x": tsp.get("x"),
                        "y": tsp.get("y"),
                    }

                # === 2c. popup ===
                if "popup" in data_record:
                    popup = data_record["popup"]
                    hmi_branch["popup"] = {
                        "name": popup.get("name"),
                        "opened": popup.get("opened"),
                    }

        # === 3. Handle mousePosition (not tied to specific flight) ===
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
        print("Flight Plan UUID: ", flight_plan_uuid)
        callsign = self.uuid_callsign_map.get(flight_plan_uuid)
        print("Callsign from UUID: ", callsign)
        
        if not callsign:
            return  # or raise a warning/log error if UUID not found

        flight_key = f"flight_{callsign}"
        predicted_points = data_record.get("points", [])

        # Ensure the structure exists and only update the "predicted" section
        route_group = self.data_repository \
            .setdefault(timestamp, {}) \
            .setdefault(flight_key, {}) \
            .setdefault("Flight", {}) \
            .setdefault("routeTrajectoryGroup", {})

        route_group["predicted"] = {
            "element": predicted_points
        } 
        
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

        Args:
            fixm_record (str): The FIXM XML string of a flight record.
            callsign (str): The aircraft ID.
        """
        self._create_new_repo_state(timestamp)

        try:
            root = ET.fromstring(fixm_record)
        except ET.ParseError as e:
            print(f"XML Parsing Error in _update_data_repository_fixm: {e}")
            return {}

        def element_to_dict(element):
            node = {}
            tag = element.tag.split("}")[-1]

            # Process attributes
            for attr_name, attr_value in element.attrib.items():
                node[attr_name] = attr_value

            # Process children elements
            children = list(element)
            if children:
                for child in children:
                    child_tag = child.tag.split("}")[-1]
                    child_dict = element_to_dict(child)

                    if child_tag not in node:
                        node[child_tag] = child_dict
                    else:
                        # Convert to list if multiple entries for the same tag
                        if not isinstance(node[child_tag], list):
                            node[child_tag] = [node[child_tag]]
                        node[child_tag].append(child_dict)
            else:
                # Leaf node: return text content directly if available
                text = element.text.strip() if element.text else ""
                if tag == "pos" and text:
                    try:
                        lat_str, lon_str = text.split()
                        return {
                            "lat": float(lat_str),
                            "lon": float(lon_str)
                        }
                    except ValueError:
                        # If the split doesn't work as expected, fallback to raw text
                        return text
                return text if text else node

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
        
        self.last_timestamp = timestamp  # Store the last processed timestamp

    def _create_new_repo_state(self, new_timestamp):
        """Create a new repository state for the next timestamp."""
        previous_state = self.data_repository.get(self.last_timestamp, {})
        self.data_repository[new_timestamp] = copy.deepcopy(previous_state)

        print(f"Created new repository state for timestamp: {new_timestamp}, used previous state: {self.last_timestamp}")
        
    def convert_event_data(self, json_record):
        """Convert non-FIXM data to RDF graph."""

        # This loop handles JSON -> RDF graph conversion TODO: refactor this into a separate function
        for key, value in json_record.items():
            print(f"Processing key: {key}, value: {value}")
            if key == "timestamp":
                timestamp_uri = URIRef(value)
                self.graph.add((timestamp_uri, RDF.type, FIXM.Timestamp))
                self._add_literal(timestamp_uri, FIXM.timestampValue, value, XSD.dateTime)

            elif key.startswith("Flight_"): # TODO: which data is this?
                self._process_flight(timestamp_uri, key, value)

            elif key == "HMI": # TODO: there is no key HMI in the original data, we should split it by alert/label_position/...
                self._process_hmi(timestamp_uri, value)

            elif key == "trajectory":
                if "flight_plan_uuid" in value.keys():
                    # TODO: join the flight_plan_uuid with the flights from the repository
                    flight_uri = ... # URIRef(FIXM[f"flight_{value['flight_plan_uuid']}"])
                    self._process_trajectory(flight_uri, value)

            elif key == "track":
                pass # TODO: implement track processing

        # This handles the data repository creation
        # TODO: can we convert repository to turtle, without additional rdf graph creation?
        self._update_data_repository_events(json_record)

        print(json.dumps(self.data_repository, indent=4))
        # print(json.dumps(self.track_number_callsign_map, indent=4))
        # Test print
        """for timestamp, data in self.data_repository.items():
            print(f"Timestamp: {timestamp}")
            for key, _ in data.items():
                if "flight" in key.lower():
                    print(f"{key}")"""
        
    def convert_fixm_data(self, xml_record):
        """Convert FIXM/XML data to RDF graph."""
        
        timestamp = xml_record["timestamp"]

        flight_plan = xml_record["flight_plan"]
        fixm_data = flight_plan["fixm"]
        uuid = flight_plan["uuid"]
        self._match_uuids_and_callsigns(uuid, fixm_data)

        # XML -> RDF graph conversion
        def process_element(element, parent_subject=None):
            tag = element.tag.split("}")[-1]
            subject_tag = tag[0].upper() + tag[1:]
            subject_id = self._generate_index(subject_tag, timestamp)

            if "base" in element.tag:
                subject = FB[subject_id]
                self.graph.add((subject, RDF.type, URIRef(FB[subject_tag])))
            else:
                subject = FIXM[subject_id]
                self.graph.add((subject, RDF.type, URIRef(FIXM[subject_tag])))

            predicates = {}
            for attr_name, attr_value in element.attrib.items():
                attr_uri = FIXM[attr_name]
                predicates[attr_uri] = [{"value": attr_value, "type": "literal"}]

            for child in element:
                child_tag = child.tag.split("}")[-1]
                if "base" in child.tag:
                    predicate_uri = FB[child_tag]
                else:
                    predicate_uri = FIXM[child_tag]

                if list(child):  # Has sub-elements
                    child_subject = process_element(child, subject)

                    if predicate_uri not in predicates:
                        predicates[predicate_uri] = [
                            {"value": child_subject, "type": "uri"}
                        ]
                elif child.text:
                    for attr_name, attr_value in child.attrib.items():
                        if "base" in child.tag:
                            attr_uri = FB[attr_name]
                        else:
                            attr_uri = FIXM[attr_name]
                        predicates[attr_uri] = [
                            {"value": attr_value, "type": "literal"}
                        ]

                    predicates[predicate_uri] = [
                        {"value": child.text.strip(), "type": "literal"}
                    ]

            if parent_subject:
                if "base" in element.tag:
                    parent_predicate = FB[tag]
                else:
                    parent_predicate = FIXM[tag]
                self.graph.add((parent_subject, parent_predicate, subject))

            self._add_triples(subject, predicates)

            return subject
        
        try:
            root = ET.fromstring(fixm_data)
        except ET.ParseError as e:
            print(f"XML Parsing Error in _update_data_repository_fixm: {e}")
            return {}
        
        process_element(root)
        print("RDF conversion completed.")

        # This handles the data repository creation
        # TODO: can we convert repository to turtle, without additional rdf graph creation?
        matched_callsign = self.uuid_callsign_map.get(uuid)
        self._update_data_repository_fixm(fixm_data, timestamp, matched_callsign)

        print(json.dumps(self.data_repository, indent=4))
        # print(json.dumps(self.track_number_callsign_map, indent=4))
        # print("Timestamp repository: ", self.data_repository)
        # print("Timestamp track_num to callsign map: ", self.track_number_callsign_map)
        # Test print
        """for timestamp, data in self.data_repository.items():
            print(f"Timestamp: {timestamp}")
            for key, _ in data.items():
                if "flight" in key.lower():
                    print(f"{key}")"""

    def serialize(self, format="turtle"):
        """Serialize the RDF graph to a string/turtle format."""
        return self.graph.serialize(format=format)

    # For testing purposes, you can save the RDF graph to a file
    def save(self, file_path, format="turtle"):
        """Save the RDF graph to a turtle file."""
        with open(file_path, "w", encoding="utf-8") as out_file:
            out_file.write(self.serialize(format=format))