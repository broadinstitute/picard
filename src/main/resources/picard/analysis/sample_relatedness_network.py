
import pandas as pd
import numpy as np
import argparse
from graphviz import Graph
from colour import Color
import umap.umap_ as umap
import matplotlib.pyplot as plt
import matplotlib


'''
Creates visualization of sample relatedness (using LOD scores from CrosscheckFingerprints)
Samples with low coverage at fingerprinting sites will be more opaque than other samples

Run command:

python3 sample_relatedness_network.py \
	-M /path/to/matrix_output.txt \
	-S /path/to/sample_individual_map.tsv

'''

hex_codes = ["#000000","#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6", "#FFDBE5", "#0000A6","#63FFAC","#8FB0FF", "#5A0007", "#FEFFE6", "#4FC601", "#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF", "#A1C299", "#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72","#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28","#C8D0F6","#A3A489","#806C66","#222800","#BF5650","#E83000","#66796D","#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51","#C895C5","#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94","#7ED379","#012C58","#7A7BFF","#D68E01","#353339","#78AFA1","#FEB2C6","#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5","#9556BD","#6A714A","#001325","#02525F","#0AA3F7","#E98176","#DBD5DD","#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75","#8D8546","#9695C5","#E773CE","#D86A78","#3E89BE","#CA834E","#518A87","#5B113C","#55813B","#E704C4","#00005F","#A97399","#4B8160","#59738A","#FF5DA7","#F7C9BF","#643127","#513A01","#6B94AA","#51A058","#A45B02","#1D1702","#E20027","#E7AB63","#4C6001","#9C6966","#64547B","#97979E","#006A66","#391406","#F4D749","#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4","#00D891","#15A08A","#BC65E9","#FFFFFE","#C6DC99","#203B3C","#671190","#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527","#8BB400","#797868","#C6005A","#3B000A","#C86240","#29607C","#402334","#7D5A44","#CCB87C","#B88183","#AA5199","#B5D6C3","#A38469","#9F94F0","#A74571","#B894A6","#71BB8C","#00B433","#789EC9","#6D80BA","#953F00","#5EFF03","#E4FFFC","#1BE177","#BCB1E5","#76912F","#003109","#0060CD","#D20096","#895563","#29201D","#5B3213","#A76F42","#89412E","#1A3A2A","#494B5A","#A88C85","#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F","#BDC9D2","#9FA064","#BE4700","#658188","#83A485","#453C23","#47675D","#3A3F00","#061203","#DFFB71","#868E7E","#98D058","#6C8F7D","#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46","#D25B88","#5B656C","#00B57F","#545C46","#866097","#365D25","#252F99","#00CCFF","#674E60","#FC009C","#92896B"]

def extract_sample_id(file_string):
	#extract sample name from file name
	sample_string = file_string.rsplit(':', 1)[-1]
	return sample_string

def convert_LOD_to_float(LOD):
	LOD = str(LOD)
	LOD_no_comma = LOD.replace(',', '')
	LOD_float = float(LOD_no_comma)
	return LOD_float

def LOD_to_distance(LOD):
	# adjust parameters here to change graph structure
	if (LOD > -10):
		distance = 1 / (1 + np.exp(-.5 * (-LOD + 4)))
	else:
		distance = 1.2
	return distance

def reformat_matrix_output(matrix_output):
	# Change file names to sample names
	# Also change LOD scores to distances for graphing 
	matrix_df = pd.read_csv(matrix_output, sep = "\t", error_bad_lines=False)
	matrix_df['FILE'] = matrix_df['FILE'].apply(extract_sample_id)
	new_names = {}
	for column in matrix_df:
		new_names[column] = extract_sample_id(column)
		if column != 'FILE':
			matrix_df[column] = matrix_df[column].apply(convert_LOD_to_float)
	matrix_df = matrix_df.rename(columns=new_names)
	return matrix_df

def rescale_matrix(matrix_df):
	new_matrix_df = pd.DataFrame()
	for column in matrix_df:
	    if column == 'FILE':
	        continue
	    new_matrix_df[column] = matrix_df[column].apply(LOD_to_distance)
	return new_matrix_df

def create_sample_map(sample_individual_map_df):
	#Creates dictionary mapping sample name to individual ID for easy access later on
	sample_map_dict = {}
	for i, row in sample_individual_map_df.iterrows():
		sample_map_dict[row.sample_id] = row.individual
	return sample_map_dict

def create_color_map(sample_ind_map_dict, matrix_df):
	# create list of hex codes with one for each sample. Samples that were supposedly 
	# from the same individual will have the same color
	individuals = []
	for sample, individual in sample_ind_map_dict.items():
		if individual in individuals:
			continue
		individuals.append(individual)
	color_map = {}
	for i in range(0, len(individuals)):
		color_map[individuals[i]] = hex_codes[i]
	colors = []
	alphas = []
	for column in matrix_df:
		if column == 'FILE':
			continue
		mag_column = matrix_df[column].abs()
		if ((mag_column.mean()) >= 5):
			alphas.append(.5)
		else:
			alphas.append(.99)
		individual = sample_ind_map_dict[column]
		colors.append(color_map[individual])
	alphas = list(map(float, alphas))
	i = 0
	color_array = []
	for c, a in zip(colors,alphas):
		col = list(matplotlib.colors.to_rgb(c))
		col.append(a)
		color_array.append(col)
	color_array = np.array(color_array)
	return color_array

def font_color(hex_code):
	c = Color(hex_code)
	if c.luminance <= .4:
		return "white"
	else:
		return "black"


def draw_graph(matrix_df, color_array):
	reducer = umap.UMAP(n_components=2,
                    n_neighbors=2,
                    spread=3,
                    min_dist=2,
                    metric='precomputed',
                    disconnection_distance=20,
                    random_state=0)
	embeddings = reducer.fit_transform(matrix_df)
	plt.scatter(
	    embeddings[:, 0],
	    embeddings[:, 1],
	    s = 100,
	    c= color_array)
	plt.title('UMAP projection of Clusters', fontsize=24)
	plt.show()

def graph_network(sample_individual_map, matrix_output):
	sample_individual_map_df = pd.read_csv(sample_individual_map, sep = "\t", error_bad_lines=False)
	sample_ind_map_dict = create_sample_map(sample_individual_map_df)
	matrix_df = reformat_matrix_output(matrix_output)
	print(LOD_to_distance(0))
	color_array = create_color_map(sample_ind_map_dict, matrix_df)
	matrix_df = rescale_matrix(matrix_df)
	draw_graph(matrix_df, color_array)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-S', default=None, type=str, help='sample individual map')
	parser.add_argument('-M', default=None, type=str, help='matrix output from crosscheck')
	args = parser.parse_args()
	graph_network(args.S, args.M)

