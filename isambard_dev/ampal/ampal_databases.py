import json
import os

from sqlalchemy import create_engine, Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from settings import global_settings


ce_path = os.path.join(global_settings['package_path'], 'ampal', 'chemical_elements.json')
with open(ce_path, 'r') as inf:
    element_data = json.loads(inf.read())

try:
    ampal_data_engine = create_engine('sqlite:///' + os.path.join(
            global_settings['package_path'], 'ampal', 'ampal_data.db'), echo=False)
    AmpalDataDecBase = declarative_base()
    AmpalDataSession = sessionmaker(bind=ampal_data_engine)
    ampal_data_session = AmpalDataSession()
finally:
    ampal_data_session.close()


class PDBColFormat(AmpalDataDecBase):
    __tablename__ = 'pdb_col_format'

    id = Column(Integer, primary_key=True)
    atom_name = Column(String)
    atom_col = Column(String)

    def __repr__(self):
        return "<Column Format(Name='{}', Column='{}')>".format(self.atom_name, self.atom_col)


global_settings['ampal'] = {'pdb_atom_col_dict': {
        e.atom_name: e.atom_col for e in ampal_data_session.query(PDBColFormat).all()}
    }


__author__ = "Christopher W. Wood"
