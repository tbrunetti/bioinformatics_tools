import json
import csv
import sys
from pprint import pprint

def jsonToCSV():
	#note, not json.loads(); -->used for json formatted string, not json file object
	with open(sys.argv[1]) as jsonData:
		data=json.load(jsonData)
		pprint(data)
	
	#for key in data:
	#	print key

if __name__=='__main__':
	jsonToCSV();