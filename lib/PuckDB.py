import sqlite3
import os
import sys

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')


class puck_db(object):

    def __init__(self, data_base):

        self.data_base = data_base

    def get_all_pucks_from_DB(self):
        print('hello')

    def insert_mounted_crystals(self):
        print('hello')

    def kill_mounted_cryatsl(self):
        print('hello')

    def change_puck_location(self):
        print('hello')

    def execute_sqlite(self):
        print('hello')
