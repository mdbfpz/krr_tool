from rdflib import RDF, Graph, Literal, Namespace, URIRef
from rdflib.namespace import RDFS, XSD

EX = Namespace("http://example.org/")
DCTERMS = Namespace("http://purl.org/dc/terms/")

class RDFConverter:
    def __init__(self):
        """Initialize the RDFConverter."""

        self._initialize_graph()

    def _initialize_graph(self):
        """Initialize the RDF graph with basic namespaces."""

        self.graph = Graph()  # Reinitialize the graph
        self.graph.bind("rdf", RDF)
        self.graph.bind("rdfs", RDFS)
        self.graph.bind("xsd", XSD)
        self.graph.bind("ex", EX)
        self.graph.bind("dcterms", DCTERMS)

    def _add_triples(self, subject, predicates):
        """Add triples to the RDF graph."""

        for predicate_uri, objects in predicates.items():
            predicate = URIRef(predicate_uri)
            for obj in objects:
                if obj["type"] == "literal":
                    if "lang" in obj:  # If the language is specified, no datatype should be provided
                        literal = Literal(obj["value"], lang=obj["lang"])
                    else:  # If no language is specified, a datatype can be provided
                        literal = Literal(obj["value"], datatype=XSD.string)
                    self.graph.add((subject, predicate, literal))
                else:
                    obj_uri = URIRef(obj["value"])
                    self.graph.add((subject, predicate, obj_uri))

    def convert_to_rdf(self, json_data):
        """Convert the JSON data to RDF format."""

        self._initialize_graph()
        self.json_data = json_data

        for subject_uri, predicates in self.json_data.items():
            subject = URIRef(subject_uri)
            self._add_triples(subject, predicates)
        
        return self.graph

    def print_graph(self):
        """Print the RDF graph in Turtle format."""
        
        print(self.graph.serialize(format="turtle"))


    def run(self, json_data):
        """ --- """

        # Sample JSON data
        sample_json = {
            "http://example.org/about": {
                "http://purl.org/dc/terms/title": [
                    {"value": "Anna's Homepage", "type": "literal", "lang": "en"}
                ]
            }
        }

        # Convert JSON to RDF and save to file
        self.convert_to_rdf(sample_json) # TODO: pass 'json_data' as an argument here once data fetching is implemented
        print("RDF conversion completed.")

        return self.graph