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
            'xs': XS,
            'fb': FIXM
        }
    
    def _initialize_graph(self):
        """Initialize the RDF graph with basic namespaces."""

        self.graph = Graph()
        self.graph.bind("rdf", RDF)
        self.graph.bind("rdfs", RDFS)
        self.graph.bind("xsd", XSD)
        self.graph.bind("fixm", FIXM)  # Bind FIXM namespace
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

    def _to_turtle(self):
        """Convert the RDF graph into the Turtle format."""
        
        self.graph = self.graph.serialize(format="turtle")

    def _xml_to_rdf(self, xml_data):
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

        # Debug: Print the XML to be parsed
        print("Parsing XML data:\n", xml_data)

        # Parse the XML data, using the namespaces defined in the class
        parser = ET.XMLParser(encoding="utf-8")
        root = ET.fromstring(xml_data, parser=parser)

        # Process the root element and its children
        for element in root.findall(".//xs:complexType", self.namespaces):  # Use the defined namespaces
            # Generate a subject based on the element's tag and namespace
            element_tag = element.attrib.get("name", "UnnamedElement")  # Use 'name' attribute for identification
            subject = URIRef(f"http://www.fixm.aero/{element_tag}")

            predicates = {}

            # Handle attributes for the element
            for attr_name, attr_value in element.attrib.items():
                attr_uri = f"http://www.w3.org/2001/XMLSchema#{attr_name}"  # Example: mapping to XSD attributes
                predicates[attr_uri] = [{"value": attr_value, "type": "literal"}]

            # Process child elements within the xs:sequence
            sequence = element.find("xs:sequence", self.namespaces)
            if sequence:
                for sub_elem in sequence:
                    sub_elem_tag = sub_elem.attrib.get("name", sub_elem.tag.split('}')[-1])
                    predicate_uri = f"http://fixm.aero/{sub_elem_tag}"

                    # Handle literal or URI values based on the element's content
                    if sub_elem.text:
                        predicates[predicate_uri] = [{"value": sub_elem.text, "type": "literal"}]
                    else:
                        predicates[predicate_uri] = [{"value": f"http://fixm.aero/{sub_elem_tag}", "type": "uri"}]

            # Add the triples for the element and its children
            self._add_triples(subject, predicates)

        print("RDF conversion completed.")

    def run(self, xml_data):
        """Run the XML to RDF conversion."""

        self._xml_to_rdf(xml_data)
        self._to_turtle()
        
        return self.graph