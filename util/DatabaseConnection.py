#Author Graham Rockwell
#Church Lab

import MySQLdb

class DatabaseConnection:
    
    def __init__(self):
        self.user = ''
        self.password = ''
        self.host = ''
        self.port = ''
        self.database = ''
        
    def config(self,config):
        self.user = config["user"]
        self.password = config["password"]
        self.host = config["host"]
        self.port = int(config["port"])
        self.database = config["database"]
        
    def getConnection(self):
        connection = MySQLdb.connect(user=self.user,
                                          passwd=self.password,
                                          host=self.host,
                                          port=self.port,
                                          db=self.database)
        return connection
    
    
    