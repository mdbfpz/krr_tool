from rdflib import Graph, Literal, Namespace, URIRef
from rdflib.namespace import RDF, RDFS, XSD

EX = Namespace("http://example.org/")

class ReverseRDFConverter:
    def __init__(self):
        """Initialize the ReverseRDFConverter."""

        self.json_data = {}
        self._initialize_graph()

    def _initialize_graph(self):
        """Initialize the RDF graph with basic namespaces."""

        self.graph = Graph()  # Reinitialize the graph
        self.graph.bind("rdf", RDF)
        self.graph.bind("rdfs", RDFS)
        self.graph.bind("xsd", XSD)
        self.graph.bind("ex", EX)

    def convert_rdf_to_json(self, rdf_data):
        """Convert RDF data to JSON format."""
        
        if isinstance(rdf_data, Graph):
            self.graph = rdf_data
        else:
            self._initialize_graph()
            self.graph.parse(data=rdf_data, format='turtle')  # Assume the input is in Turtle format

        json_data = {}
        
        # Loop over each subject in the graph
        for subject in set(self.graph.subjects()):
            subject_str = str(subject)
            json_data[subject_str] = {}

            # Loop over predicates and objects for each subject
            for predicate, obj in self.graph.predicate_objects(subject):
                predicate_str = str(predicate)
                
                # Initialize the predicate entry if it doesn't exist yet
                if predicate_str not in json_data[subject_str]:
                    json_data[subject_str][predicate_str] = []

                # Determine the type of the object (URIRef or Literal)
                if isinstance(obj, Literal):
                    # Handle literals: include value, type (literal), and lang if available
                    literal_entry = {
                        "value": obj.toPython(),
                        "type": "literal"
                    }
                    if obj.language:
                        literal_entry["lang"] = obj.language

                    json_data[subject_str][predicate_str].append(literal_entry)

                elif isinstance(obj, URIRef):
                    # Handle URIs
                    json_data[subject_str][predicate_str].append({
                        "value": str(obj),
                        "type": "uri"
                    })

        self.json_data = json_data
        print("RDF to JSON conversion completed.")
    
    def run(self, rdf_graph_data):
        """ --- """

        # Convert RDF to JSON
        self.convert_rdf_to_json(rdf_graph_data)

        return self.json_data
