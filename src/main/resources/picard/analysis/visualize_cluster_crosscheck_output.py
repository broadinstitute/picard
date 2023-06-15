import pandas as pd
import argparse
from graphviz import Graph
from colour import Color



### Visualize ClusterCrosscheckMetrics

# Example command:
# python3 reformat_cluster_crosscheck_output.py 
#	--clustered_crosscheck_metrics clustered_metrics.tsv 
#	--sample_individual_map sample_individual_map.tsv 
#	--matrix matrix_output.txt 
#	--table samples.tsv 
#	--min_edge_weight 1


hex_codes = ["#000000","#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6", "#FFDBE5", "#0000A6","#63FFAC","#8FB0FF", "#5A0007", "#FEFFE6", "#4FC601", "#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF", "#A1C299", "#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72","#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28","#C8D0F6","#A3A489","#806C66","#222800","#BF5650","#E83000","#66796D","#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51","#C895C5","#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94","#7ED379","#012C58","#7A7BFF","#D68E01","#353339","#78AFA1","#FEB2C6","#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5","#9556BD","#6A714A","#001325","#02525F","#0AA3F7","#E98176","#DBD5DD","#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75","#8D8546","#9695C5","#E773CE","#D86A78","#3E89BE","#CA834E","#518A87","#5B113C","#55813B","#E704C4","#00005F","#A97399","#4B8160","#59738A","#FF5DA7","#F7C9BF","#643127","#513A01","#6B94AA","#51A058","#A45B02","#1D1702","#E20027","#E7AB63","#4C6001","#9C6966","#64547B","#97979E","#006A66","#391406","#F4D749","#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4","#00D891","#15A08A","#BC65E9","#FFFFFE","#C6DC99","#203B3C","#671190","#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527","#8BB400","#797868","#C6005A","#3B000A","#C86240","#29607C","#402334","#7D5A44","#CCB87C","#B88183","#AA5199","#B5D6C3","#A38469","#9F94F0","#A74571","#B894A6","#71BB8C","#00B433","#789EC9","#6D80BA","#953F00","#5EFF03","#E4FFFC","#1BE177","#BCB1E5","#76912F","#003109","#0060CD","#D20096","#895563","#29201D","#5B3213","#A76F42","#89412E","#1A3A2A","#494B5A","#A88C85","#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F","#BDC9D2","#9FA064","#BE4700","#658188","#83A485","#453C23","#47675D","#3A3F00","#061203","#DFFB71","#868E7E","#98D058","#6C8F7D","#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46","#D25B88","#5B656C","#00B57F","#545C46","#866097","#365D25","#252F99","#00CCFF","#674E60","#FC009C","#92896B"]

def remove_header(file):
	# Removes header from cluster crosscheck metrics file so pandas can read the TSV
	metrics = open(file, "r")
	new_file = open(f'{file}_no_header', "w")
	for i, line in enumerate(metrics):
	    if not line.startswith('#'):
	        new_file.write(line)

def extract_cluster_info(clusters_df):
	# Creates a dictionary mapping cluster names to a list of cluster members
	# also records which samples are in clusters with > 1 member (not singletons)
	cluster_dict = {}
	in_clusters = []
	for i, row in clusters_df.iterrows():
		if row.CLUSTER not in cluster_dict.keys():
			cluster_dict[row.CLUSTER] = []
		if row.CLUSTER_SIZE > 1:
			if row.RIGHT_SAMPLE not in in_clusters:
				in_clusters.append(row.RIGHT_SAMPLE)
			if row.LEFT_SAMPLE not in in_clusters:
				in_clusters.append(row.LEFT_SAMPLE)
		if row.RIGHT_SAMPLE not in cluster_dict[row.CLUSTER]:
			cluster_dict[row.CLUSTER].append(row.RIGHT_SAMPLE)
		if row.LEFT_SAMPLE not in cluster_dict[row.CLUSTER]:
			cluster_dict[row.CLUSTER].append(row.LEFT_SAMPLE)
	return cluster_dict, in_clusters

def add_singletons(in_clusters, table):
	# Finds samples that were included in crosscheck but not included in clusters (for some reason? not sure why yet) or were in a singleton cluster
	table_df = pd.read_csv(table, sep = "\t", error_bad_lines=False)
	singletons = []
	for index, row in table_df.iterrows():
		if row['entity:participant_id'] not in in_clusters:
			singletons.append(row['entity:participant_id']) 
	return singletons

def extract_sample_id(file_string):
	#extract sample name from file name
	sample_string = file_string.rsplit(':', 1)[-1]
	return sample_string

def reformat_matrix_output(matrix_output):
	#change file names to sample names in matrix output
	matrix_df = pd.read_csv(matrix_output, sep = "\t", error_bad_lines=False)
	matrix_df['FILE'] = matrix_df['FILE'].apply(extract_sample_id)
	new_names = {}
	for column in matrix_df:
		new_names[column] = extract_sample_id(column)
	matrix_df = matrix_df.rename(columns=new_names)
	return matrix_df

def find_edges(singletons, matrix_df, min_edge_weight):
	# find edges for singletons (if they exist)
	# In other words, find where singletons have a positive LOD score with another sample
	# return a list of thruples with sampleid1, sampleid2, and LOD score
	# also return a list of samples with no coverage at fingerprinting sites
	edges = []
	no_coverage_samples = []
	for index, row in matrix_df.iterrows():
		#initiate no coverage variable
		no_coverage = True
		#check if row file is in singletons
		row_sample = row.FILE
		if row_sample not in singletons:
			continue
		#check for positive values in row
		row_dict = row.to_dict()
		for column_sample, LOD_score in row_dict.items():
			LOD_score = str(LOD_score)
			if column_sample == 'FILE':
				continue
			LOD_score = LOD_score.replace(',', '')
			# if there is a non-zero LOD score, there is coverage at fingerprinting site
			if float(LOD_score) != 0:
				no_coverage = False
			# create edge if LOD score is above threshold
			if float(LOD_score) > min_edge_weight:
				edge = (row_sample, column_sample, LOD_score)
			#make sure this thruple (but w files reversed) isn't in alread in output
				if (column_sample, row_sample, LOD_score) in edges:
					continue
			#remove edges that connect a node to itself
				if column_sample == row_sample:
					continue
			#add thruple to output
				edges.append(edge)
		if no_coverage:
			no_coverage_samples.append(row_sample)
	return edges, no_coverage_samples

def create_sample_map(sample_individual_map_df):
	#Creates dictionary mapping sample name to individual ID for easy access later on
	sample_map_dict = {}
	for i, row in sample_individual_map_df.iterrows():
		sample_map_dict[row.sample_id] = row.individual
	return sample_map_dict

def create_color_map(sample_ind_map_dict):
	# Maps each individual to a color
	individuals = []
	for sample, individual in sample_ind_map_dict.items():
		if individual in individuals:
			continue
		individuals.append(individual)
	color_map = {}
	for i in range(0, len(individuals)):
		color_map[individuals[i]] = hex_codes[i]
	return(color_map)

def font_color(hex_code):
	# Assigns font color (black or white) based on luminance of background color
	c = Color(hex_code)
	if c.luminance <= .4:
		return "white"
	else:
		return "black"

def render_clusters(cluster_dict, sample_ind_map_dict, ind_color_map, 
	singletons, edges, no_coverage_samples):
	nodes_included = []
	g = Graph('G', filename='cluster.gv', 
		node_attr={'shape': 'record'}, 
		graph_attr={'pack': 'true', 'layout': 'fdp'})
	# Draw legend
	with g.subgraph(name='cluster_Legend') as c:
		for ind, color in ind_color_map.items():
			c.node(ind, 
    			shape="rectangle", 
				style = "filled", 
				color = color,
				fontcolor = font_color(color))
			c.attr(label='Legend', pack='true', bgcolor='gainsboro')
	# Draw each cluster
	for cluster, samples in cluster_dict.items():
		with g.subgraph(name=f'cluster_{cluster}') as c:
		    for sample in samples:
		    	nodes_included.append(sample)
		    	ind = sample_ind_map_dict[sample]
		    	color = ind_color_map[ind]
		    	c.node(f'{sample}', 
		    		shape="rectangle", 
		    		style = "filled", 
		    		color = color,
		    		fontcolor = font_color(color))
	    		c.attr(label=f'Cluster {cluster}')
	# Draw no-coverage cluster
	with g.subgraph(name='cluster_no_coverage') as c:
		for sample in no_coverage_samples:
			nodes_included.append(sample)
			ind = sample_ind_map_dict[sample]
			color = ind_color_map[ind]
			c.node(f'{sample}', 
				shape="rectangle", 
				style = "filled", 
				color = color,
				fontcolor = font_color(color))
		c.attr(label="No coverage at fingerprinting sites", bgcolor="red")
	# Draw singletons
	for sample in singletons:
		if sample in nodes_included: 
			continue
		ind = sample_ind_map_dict[sample]
		color = ind_color_map[ind]
		g.node(f'{sample}', 
			shape="rectangle", 
			style = "filled", 
			color = color,
			fontcolor = font_color(color))
	# Add edges
	for edge in edges:
		g.edge(edge[0], edge[1], label=edge[2])

	g.render('crosscheck_clusters', view=True)

def reformat(cluster_file, sample_individual_map, matrix_output, table, min_edge_weight):
	remove_header(cluster_file)
	new_file = f'{cluster_file}_no_header'
	clusters_df = pd.read_csv(new_file, sep = "\t", error_bad_lines=False)
	sample_individual_map_df = pd.read_csv(sample_individual_map, sep = "\t", error_bad_lines=False)
	cluster_dict, in_clusters = extract_cluster_info(clusters_df)
	singletons = add_singletons(in_clusters, table)
	matrix_df = reformat_matrix_output(matrix_output)
	edges, no_coverage_samples = find_edges(singletons, matrix_df, min_edge_weight)
	sample_ind_map_dict = create_sample_map(sample_individual_map_df)
	ind_color_map = create_color_map(sample_ind_map_dict)
	render_clusters(cluster_dict, sample_ind_map_dict, ind_color_map, singletons, edges, no_coverage_samples)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--clustered_crosscheck_metrics', default=None, type=str, help='path to cluster crosscheck metrics file output')
	parser.add_argument('--sample_individual_map', default=None, type=str, help='sample individual map')
	parser.add_argument('--matrix', default=None, type=str, help='matrix output from crosscheck')
	parser.add_argument('--table', default=None, type=str, help='table with file URLs to map to sample names (if necessary)')
	parser.add_argument('--min_edge_weight', default=None, type=int, help='minimum edge weight to display for singletons')
	args = parser.parse_args()
	reformat(args.clustered_crosscheck_metrics, args.sample_individual_map, args.matrix, args.table, args.min_edge_weight)
