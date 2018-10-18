import os, sys

#runs SpoTyping script on all paired fastq files
def main(argv):
	forward = []
	reverse = []
	for file in os.listdir("../../Meehan_TB/"):
		if(file.endswith("_1.fastq.gz")):
			forward.append(file)
		elif(file.endswith("_2.fastq.gz")): 
			reverse.append(file)
	forward.sort()
	reverse.sort()
	print(forward[0], reverse[0])
	
	i = len(forward)
	for j in range(i):
		os.system('python SpoTyping.py ../../Meehan_TB/%s ../../Meehan_TB/%s -o spotyping_output/spo.out' % (forward[j],reverse[j]))

	
if __name__ == "__main__":
   main(sys.argv[1:])
