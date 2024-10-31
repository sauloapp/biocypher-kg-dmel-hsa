# Author Abdulrahman S. Omar <xabush@singularitynet.io>

# class Adapter:
#     def __init__(self, write_properties, add_provenance):
#         self.write_properties = write_properties
#         self.add_provenance = add_provenance
#     def get_nodes(self):
#         pass

#     def get_edges(self):
#         pass


from abc import ABC, abstractmethod

# Adapter now inherits from ABC, making it an abstract class
class Adapter(ABC):
    def __init__(self, write_properties, add_provenance):
        """
        Constructor for the Adapter class.
        """
        self.write_properties = write_properties
        self.add_provenance = add_provenance

    @abstractmethod
    def get_nodes(self):
        """Abstract method to be implemented by subclasses."""
        pass

    @abstractmethod
    def get_edges(self):
        """Abstract method to be implemented by subclasses."""
        pass


# EdgeAdapter now inherits from Adapter
class EdgeAdapter(Adapter):
    
    def __init__(self, write_properties, add_provenance, table_file_name):
        """
        Constructor for the EdgeAdapter class.
        Calls the Adapter constructor and initializes additional properties.
        """
        super().__init__(write_properties, add_provenance)
        self.table_file_name = table_file_name

    def get_edges(self):
        """
        Concrete method that subclasses can use or extend.
        Iterates over rows in the data table and generates edges.
        """
        self.data_table = self.get_table(self.table_file_name)
        for row in self.data_table:
            edge_source_id = self.get_source_id(row)
            edge_target_id = self.get_target_id(row)
            if self.valid_ids(edge_source_id, edge_target_id):
                properties = self.get_properties(row)
                if self.valid_properties(properties):
                    yield edge_source_id, edge_target_id, self.label, properties

    @abstractmethod
    def get_table(self, table_file_name):
        """Abstract method that must be implemented by subclasses."""
        pass

    @abstractmethod
    def get_source_id(self, row):
        """Abstract method that must be implemented by subclasses."""
        pass

    @abstractmethod
    def get_target_id(self, row):
        """Abstract method that must be implemented by subclasses."""
        pass

    @abstractmethod
    def valid_ids(self, source_id, target_id):
        """Abstract method that must be implemented by subclasses."""
        pass

    @abstractmethod
    def get_properties(self, row):
        """Abstract method that must be implemented by subclasses."""
        pass

    @abstractmethod
    def valid_properties(self, properties):
        """Abstract method that must be implemented by subclasses."""
        pass




# EdgeAdapter now inherits from Adapter
class NodeAdapter(Adapter):