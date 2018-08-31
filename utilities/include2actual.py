import os
import sys

def get_includes(filename,inclst):
	
	list = []
	
	directory = filename.split("/")[:-1]
	
	dirtext = ""
						
	for d in directory:
		dirtext = dirtext + d + "/"
	
	with open(filename,"r") as f:
		for line in f:

			if not line.startswith("//") and 'include "' in line and '#' in line:
				newfile = line.split('"')[1]
			
				if os.path.exists(dirtext + newfile):
					if dirtext + newfile not in inclst:
						list.append(get_includes(dirtext + newfile,inclst))
						inclst.append(dirtext + newfile)

				else:
					list.append(line)
			else:
				list.append(line)
	return list

def write_file(fhandle,list):

	for row in list:
		if type(row) is str:
			fhandle.write(row);
		
		else:
			write_file(fhandle,row)
	fhandle.write("\n")
		
if __name__ == '__main__':

	if len(sys.argv) > 1:
		result = get_includes(sys.argv[1],[])
		
		with open("combine.cpp","w") as f:
			write_file(f,result)
			pass
		