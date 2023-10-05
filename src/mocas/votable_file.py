from astropy.io.votable.tree import VOTableFile, Resource, Table as VOTable, Field


class TAPUploadVOTableFile(VOTableFile):

    def __init__(self, rows: list, field_definitions: dict, **kwargs):
        """
        Create a VOTable given a list of rows and a dictionary of field definitions for those rows.

        :param rows
        :param field_definitions: a dictionary of column definitions
        :return:
        """
        super().__init__(**kwargs)

        # ...with one resource...
        resource = Resource()
        self.resources.append(resource)

        # ... with one table
        table = VOTable(self)
        resource.tables.append(table)

        # Define some fields
        table.fields.extend([
            Field(self, name=name, **field_definitions[name]) for name in field_definitions])

        table.create_arrays(len(rows))
        for i in range(len(rows)):
            table.array[i] = rows[i]


