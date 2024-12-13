import xml.etree.ElementTree as ET
from rdflib import RDF, Graph, Literal, Namespace, URIRef
from rdflib.namespace import RDFS, XSD

FIXM = Namespace("http://www.fixm.aero/fixm/4.2.0")
XS = Namespace("http://www.w3.org/2001/XMLSchema")  # Add the xs namespace

class RDFConverter:
    def __init__(self):
        """Initialize the RDFConverter."""

        self._initialize_graph()
        # Define the namespaces used in the XML data
        self.namespaces = {
            "xs": "http://www.w3.org/2001/XMLSchema",
            "fb": "http://www.fixm.aero/fixm/4.2.0"
        }
        self.initial_flight_plans = {} # For each flight key, there is a list of points = flight plan
        self.updated_flight_plans = {} # For each flight key, there is an updated list of points = flight plan
    
    def _extract_flight_plans(self, xml_data):
        """
        For each flight in received flight plans in FIXM/XML format, save its flight plan as a
        list of points on the route.
        """

        plans = {}

        for flight_plan in xml_data:
            plans[flight_plan] = ...
                
        return plans
    
    def _update_flight_plans(self, extracted_plans):
        """
        Updates the flight plans by setting initial flight plans if they do not exist, 
        or updated flight plans if initial plans are already present.
        """

        if not self.initial_flight_plans:
            self.initial_flight_plans = extracted_plans
        else:
            self.updated_flight_plans = extracted_plans
    
    def _initialize_graph(self):
        """Initialize the RDF graph with basic namespaces."""

        self.graph = Graph()
        self.graph.bind("rdf", RDF)
        self.graph.bind("rdfs", RDFS)
        self.graph.bind("xsd", XSD)
        self.graph.bind("fixm", FIXM)
        self.graph.bind("xs", XS)

    def _add_triples(self, subject, predicates):
        """Add triples to the RDF graph."""

        for predicate_uri, objects in predicates.items():
            predicate = URIRef(predicate_uri)
            for obj in objects:
                if obj["type"] == "literal":
                    if "lang" in obj:
                        literal = Literal(obj["value"], lang=obj["lang"])
                    else:
                        literal = Literal(obj["value"], datatype=XSD.string)
                    self.graph.add((subject, predicate, literal))
                else:
                    obj_uri = URIRef(obj["value"])
                    self.graph.add((subject, predicate, obj_uri))

    def _rdf_to_turtle(self):
        """Convert the RDF graph into the Turtle format."""

        return self.graph.serialize(format="turtle")

    def _xml_to_rdf(self, xml_data: str):
        """Convert the FIXM XML data to RDF format."""

        self._initialize_graph()

        # Ensure the XML data has the correct declaration if necessary
        if not xml_data.startswith("<?xml"):
            xml_data = '<?xml version="1.0" encoding="UTF-8"?>\n' + xml_data

        # Add the namespace declarations to the XML string
        xml_data = xml_data.replace(
            "<xs:complexType",
            "<xs:complexType xmlns:xs='http://www.w3.org/2001/XMLSchema' xmlns:fb='http://www.fixm.aero/fixm/4.2.0'"
        )

        # Debugging: Print the XML to be parsed
        print("Modified XML data:\n", xml_data)

        # Parse the XML data
        parser = ET.XMLParser(encoding="utf-8")
        root = ET.fromstring(xml_data, parser=parser)

        # Debug: Print the root tag and check namespaces
        print(f"Root tag: {root.tag}")
        print(f"Namespaces used: {self.namespaces}")

        # Since root is already a complexType element, process it directly
        element_tag = root.attrib.get("name", "UnnamedElement")
        subject = URIRef(f"http://www.fixm.aero/{element_tag}")

        predicates = {}

        # Handle attributes for the root element
        for attr_name, attr_value in root.attrib.items():
            attr_uri = XS[attr_name]  # Use XS namespace for attributes
            predicates[attr_uri] = [{"value": attr_value, "type": "literal"}]

        # Process child elements within the xs:sequence
        sequence = root.find("xs:sequence", self.namespaces)
        if sequence:
            for sub_elem in sequence:
                sub_elem_tag = sub_elem.attrib.get("name", sub_elem.tag.split('}')[-1])
                predicate_uri = FIXM[sub_elem_tag]  # Use FIXM namespace for elements

                # Instead of using the element name as the object, use the element's type as the URI
                sub_elem_type = sub_elem.attrib.get("type")
                if sub_elem_type:
                    object_uri = URIRef(f"http://fixm.aero/{sub_elem_type.split(':')[-1]}")
                    predicates[predicate_uri] = [{"value": object_uri, "type": "uri"}]
                else:
                    # If no type is specified, use the element name as a literal value
                    predicates[predicate_uri] = [{"value": sub_elem_tag, "type": "literal"}]

        # Add the triples for the root element and its children
        self._add_triples(subject, predicates)

        print("RDF conversion completed.")

    def run(self, xml_data):
        """Run the XML to RDF conversion."""

        #TODO: update flight plans here

        self._xml_to_rdf(xml_data)
        print(f"Graph contains {len(self.graph)} triples after XML to RDF conversion.")
        print("Triples in the graph:")
        for s, p, o in self.graph:
            print(f"Subject: {s}, Predicate: {p}, Object: {o}")

        print("Turtle format: ", self._rdf_to_turtle())

        return self._rdf_to_turtle()  # Return the serialized graph in Turtle format