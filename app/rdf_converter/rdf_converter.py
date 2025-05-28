from rdflib import Graph, URIRef, Literal, Namespace
from rdflib.namespace import RDF, XSD
import json

# Define namespaces
FIXM = Namespace("http://www.fixm.aero/flight/4.3/")
FB = Namespace("http://www.fixm.aero/base/4.3")
XS = Namespace("http://www.w3.org/2001/XMLSchema")
HMI = Namespace("http://www.fixm.aero/hmi/1.0/") # TODO: this should not be a part of FIXM, but a separate namespace, maybe 'plain' extension

class RDFConverter:
    def __init__(self):
        """Initialize the RDF graph and bind namespaces."""
        self.graph = Graph()
        self.graph.bind("fx", FIXM)
        self.graph.bind("hmi", HMI)
        self.graph.bind("xs", XS)

    def _add_literal(self, subject, predicate, value, datatype=XSD.string):
        """Add a literal triple to the graph."""
        self.graph.add((subject, predicate, Literal(value, datatype=datatype)))

    def _add_uri(self, subject, predicate, obj):
        """Add a URI triple to the graph."""
        self.graph.add((subject, predicate, obj))

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
            lat_lon_str = " ".join(
                f"{point['lat_deg']} {point['lon_deg']}" for point in trajectory_data["points"]
            )
            self._add_literal(trajectory_uri, FIXM.trajectoryPoints, lat_lon_str)

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

    def convert(self, data):
        """Convert JSON data to RDF graph."""
        for timestamp, content in data.items():
            timestamp_uri = URIRef(timestamp)
            self.graph.add((timestamp_uri, RDF.type, FIXM.Timestamp))
            self._add_literal(timestamp_uri, FIXM.timestampValue, timestamp, XSD.dateTime)

            for entity_key, entity_data in content.items():
                if entity_key.startswith("Flight_"):
                    self._process_flight(timestamp_uri, entity_key, entity_data)
                elif entity_key == "HMI":
                    self._process_hmi(timestamp_uri, entity_data)

        return self.graph

    def serialize(self, format="turtle"):
        """Serialize the RDF graph to a string/turtle format."""
        return self.graph.serialize(format=format)

    # For testing purposes, you can save the RDF graph to a file
    def save(self, file_path, format="turtle"):
        """Save the RDF graph to a turtle file."""
        with open(file_path, "w", encoding="utf-8") as out_file:
            out_file.write(self.serialize(format=format))


"""if __name__ == "__main__":
    json_file_path = "dictModel.json"
    output_file_path = "from_dict_to_rdf.ttl"

    try:
        with open(json_file_path, 'r', encoding='utf-8') as file:
            data = json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        data = {}

    converter = RDFConverter()
    converter.convert(data)
    converter.save(output_file_path)
    print(f"RDF data has been saved to {output_file_path}")"""