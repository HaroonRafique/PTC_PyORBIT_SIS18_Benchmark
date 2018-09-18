from itertools import chain

def read_parameter_file(filename):
	file_in = open(filename,"r")
	input_lines = file_in.readlines()
	file_in.close()

	input_lines_split = map(lambda i: input_lines[i].replace('=','').replace(';','').split(), range(len(input_lines)))
	#remove commented lines
	for i in input_lines_split[::-1]:
		if '!' in ''.join(i):
			input_lines_split.remove(i)
			# print 'removed line: %s'%i
	input_list = list(chain.from_iterable(input_lines_split))
	parameters = input_list[::2]
	values = [float(i) for i in input_list[1::2]]
	input_dictionary = dict(zip(parameters, values))
	return input_dictionary