from sys import argv

def check(file_name,num):

	hand = open(file_name,'r')
	for l in hand:
		if l[0] != '#':
			alt = l.split('\t')[4]
			alts = alt.split(',')
			if len(alts) >= num :
				print l.strip()
	
	hand.close()

def main():
	num = 2
	if argv[2] != None:
		num = int(argv[2])
	
	file_name = argv[1]
	check(file_name,num)

if __name__=="__main__":
	main()
