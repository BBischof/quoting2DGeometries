import sys, json, math
from itertools import izip, islice, product
import networkx as nx 

PADDING = 0.1 #inches
MATERIAL_COST = 0.75 #dollars/in^2
LASER_SPEED = 0.5 #in/s
TIME_COST = 0.07 #dollars/s
PI = math.pi

def costFormat(x): 
	'''Coverts numbers to a number format that makes cents.'''
	return "$"+ "{:.2f}".format(float(x)) 

def createGraph(schema): 
	'''Creates a graph from the points and edges.''' 
	g = nx.Graph()
	g.add_nodes_from([int(v) for v in schema["Vertices"].keys()])
	g.add_edges_from([tuple(schema["Edges"][e]["Vertices"]) for e in schema["Edges"].keys()])
	return g

def pointLookup((x,y), schema): 
	'''Helper for reverse point lookup by coords.'''
	return [s for s in schema["Vertices"].keys() 
		if schema["Vertices"][s]["Position"]["X"] == x 
		and schema["Vertices"][s]["Position"]["Y"] == y][0]

def edgeLookup((x,y), schema):
	'''Helper for reverse edge lookup by vertices.'''
	return [s for s in schema["Edges"].keys()
		if int(x) in schema["Edges"][s]["Vertices"]
		and int(y) in schema["Edges"][s]["Vertices"]][0]

def findPath(schema_graph):
	'''Determines if the schema contains a closed curve.'''
	cycles = nx.cycle_basis(schema_graph)
	if len(cycles) > 0:
		return cycles
	else:
		raise nx.exception.NetworkXNoCycle('No closed cycle found.')

def validateSchema(name, schema): 
	'''Checks to make sure the schema is as expected, computes a graph from the data, 
	and returns the set of paths in the schema and the graph.'''
	if "Edges" not in schema.keys() or "Vertices" not in schema.keys():
		raise ValueError(name + " is not a valid schema")
	if (len(schema["Edges"].keys()) != len(set(schema["Edges"].keys())) 
		or len(schema["Vertices"].keys()) != len(set(schema["Vertices"].keys()))):
		raise ValueError(name + " is not a valid schema")
	g = createGraph(schema)
	return findPath(g), g

def computeValues(schema): 
	'''Stores the relevant values here to avoid multiple calculations.
	Organized for readability.'''
	values_dict = {}
	for edge_id in schema["Edges"].keys():
		edge = schema["Edges"][edge_id]
		start_id = ""
		start_coords = (0.0, 0.0)
		finish_id = ""
		finish_coords = (0.0, 0.0)
		dist = 0.0
		center = (0.0,0.0)
		radius = 0.0
		chord_length = 0.0
		chord_midpoint = (0.0,0.0)
		sector_length = 0.0
		sector_height = 0.0
		box_points = []
		if edge["Type"] == "LineSegment":
			start_id = str(edge["Vertices"][0])
			finish_id =  str(edge["Vertices"][1])
			start_coords = (schema["Vertices"][start_id]["Position"]["X"], schema["Vertices"][start_id]["Position"]["Y"])
			finish_coords = (schema["Vertices"][finish_id]["Position"]["X"], schema["Vertices"][finish_id]["Position"]["Y"])
			dist = math.sqrt((finish_coords[0] - start_coords[0])**2 
				+ (finish_coords[1] - start_coords[1])**2)
			values_dict[edge_id]= {"start_id": start_id, 
				"finish_id": finish_id, 
				"dist": dist}
		else:
			start_id = str(edge["ClockwiseFrom"])
			finish_id = str([p for p in edge["Vertices"] if str(p) != start_id][0])
			start_coords = (schema["Vertices"][start_id]["Position"]["X"], schema["Vertices"][start_id]["Position"]["Y"])
			finish_coords = (schema["Vertices"][finish_id]["Position"]["X"], schema["Vertices"][finish_id]["Position"]["Y"])
			center = (edge["Center"]["X"], edge["Center"]["Y"])
			radius = math.sqrt((center[0] - start_coords[0])**2 
				+ (center[1] - start_coords[1])**2)
			dx = finish_coords[0] - start_coords[0]
			dy = finish_coords[1] - start_coords[1]
			chord_length = math.sqrt(dx**2 + dy**2)
			chord_midpoint = ((start_coords[0]+finish_coords[0])/2,
				(start_coords[1]+finish_coords[1])/2)
			sector_length = 2*radius*math.asin(chord_length/(2*radius))
			sector_height = radius - math.sqrt(radius**2-chord_length**2/4)
			box_points = [(start_coords[0] + sector_height*(dy)/chord_length, start_coords[1] - sector_height*(dx)/chord_length),
				(start_coords[0] - sector_height*(dy)/chord_length, start_coords[1] + sector_height*(dx)/chord_length),
				(finish_coords[0] - sector_height*(dy)/chord_length, finish_coords[1] + sector_height*(dx)/chord_length),
				(finish_coords[0] + sector_height*(dy)/chord_length, finish_coords[1] - sector_height*(dx)/chord_length) ]
			values_dict[edge_id] = {"start_id": start_id, 
				"finish_id": finish_id, 
				"center": center, 
				"radius": radius, 
				"chord_length": chord_length, 
				"chord_midpoint": chord_midpoint, 
				"sector_length": sector_length, 
				"sector_height": sector_height, 
				"box_points": box_points}
	return values_dict


def arcLength(schema, values_dict): 
	'''Computes the price of cutting out the specified arc.'''
	price = 0.0
	for edge_id in schema["Edges"].keys():
		edge = schema["Edges"][edge_id]
		if edge["Type"] == "LineSegment":
			price += (values_dict[edge_id]["dist"]
						/LASER_SPEED*TIME_COST)
		else:
			price += (values_dict[edge_id]["sector_length"]
						/(LASER_SPEED*math.exp(-1/values_dict[edge_id]["radius"]))*TIME_COST)
	return price

def turn(p, q, r): 
	'''From https://gist.github.com/tixxit/242402.'''
	return cmp((q[0] - p[0])*(r[1] - p[1]) - (r[0] - p[0])*(q[1] - p[1]), 0)

def keep_left(hull, r): 
	'''From https://gist.github.com/tixxit/242402.'''
	while len(hull) > 1 and turn(hull[-2], hull[-1], r) != 1:
	        hull.pop()
	if not len(hull) or hull[-1] != r:
	    hull.append(r)
	return hull

def grahams_hull(points): 
	'''From https://gist.github.com/tixxit/242402.'''
	points = sorted(points)
	l = reduce(keep_left, points, [])
	u = reduce(keep_left, reversed(points), [])
	return l.extend(u[i] for i in xrange(1, len(u) - 1)) or l

def computeInnerHull(schema, schema_graph):
	'''helper function to calculate an inital hull based on the points in a cycle'''
	hull_points = []
	hulls = []
	for path in findPath(schema_graph):
		hull_points.append([pointLookup(h, schema) 
			for h in grahams_hull([(schema["Vertices"][str(p)]["Position"]["X"], schema["Vertices"][str(p)]["Position"]["Y"]) 
			for p in path])])
	return hull_points

def computeConvexHull(schema, schema_graph, values_dict):
	'''Computes the convex hull by starting with an initial hull, and then expanding
	if there are some extrusions to consider.'''
	initial_hull = computeInnerHull(schema, schema_graph)
	final_hull = []
	if len([e for e in schema["Edges"].keys() if schema["Edges"][e]["Type"] == "CircularArc"]) == 0:
		for hull in initial_hull:
			final_hull.append([(schema["Vertices"][str(h)]["Position"]["X"], schema["Vertices"][str(h)]["Position"]["Y"]) 
				for h in hull])
		return final_hull
	else: 
		'''Only goes through this if there exists a circular arc.'''
		for hull in initial_hull:
			new_hull = [(schema["Vertices"][str(h)]["Position"]["X"], schema["Vertices"][str(h)]["Position"]["Y"]) 
				for h in hull]
			points = hull+[hull[0]] #let's not talk about this hack...
			for a,b in izip(points, islice(points, 1, None)):
				path = [p for p in nx.shortest_path(schema_graph, int(a), int(b))]
				arcs_on_path = [arc for arc in [edgeLookup((x,y), schema) 
					for x,y in izip(path, islice(path, 1, None))] 
					if schema["Edges"][arc]["Type"] == "CircularArc"]
				extrusions = [x for x in arcs_on_path #Dirty trick to check if an arc is an extrusion from the hull.
					if path.index(int(values_dict[x]["start_id"])) > path.index(int(values_dict[x]["finish_id"]))]
				for x in extrusions:
					new_hull.extend(values_dict[x]["box_points"])
			final_hull.append(grahams_hull(new_hull))
		return final_hull

def findMinimumRectangle(hull):
	'''From convex hull, we compute the areas of all bounding boxes. We incorporate PADDING here.'''
	p = hull+[hull[0]]
	caliper_sizes = {}
	for s in product([x for x in product(izip(p, islice(p, 1, None)), hull) 
		if x[1] not in x[0]], product(hull,hull)):
	    x1, y1 = s[0][0][0][0], s[0][0][0][1]
	    x2, y2 = s[0][0][1][0], s[0][0][1][1]
	    x3, y3 = s[0][1][0], s[0][1][1]
	    x4, y4 = s[1][0][0], s[1][0][1]
	    x5, y5 = s[1][1][0], s[1][1][1]
	    height_numerator = ((y2-y1)*x3-(x2-x1)*y3+x2*y1-y2*x1) #avoid computing a sqrt for every single rect
	    width_numerator = ((x5-x4)*(x2-x1)+(y5-y4)*(y2-y1))
	    segment_length_square = ((x2-x1)**2+(y2-y1)**2)
	    area = math.fabs(height_numerator*width_numerator)/segment_length_square
	    if s[0][0] in caliper_sizes.keys():
	        if caliper_sizes[s[0][0]][0] < area:
	            caliper_sizes[s[0][0]] = (area, (height_numerator, width_numerator, segment_length_square))
	    else:
	        caliper_sizes[s[0][0]] = (area, (height_numerator, width_numerator, segment_length_square))
	best_edge = min(caliper_sizes, key=caliper_sizes.get)        
	return ((math.fabs(caliper_sizes[best_edge][1][0])/math.sqrt(caliper_sizes[best_edge][1][2])+PADDING)
		*(math.fabs(caliper_sizes[best_edge][1][1])/math.sqrt(caliper_sizes[best_edge][1][2])+PADDING))

def rotateCalipers(schema, schema_graph, values_dict):
	'''Creates hulls, computes areas, and optimizes.'''
	hulls = computeConvexHull(schema, schema_graph, values_dict)
	price = 0.0
	areas = []
	for hull in hulls:
		areas.append(findMinimumRectangle(hull))
	for region in areas:
		price += MATERIAL_COST*region
	return price

def main(schemas):
	values_dict = {}
	for s in schemas:
		with open(s) as data_file:  
			try:
				data = json.load(data_file)
				path, graph = validateSchema(s, data)
				values = computeValues(data)
				print s, (costFormat(arcLength(data, values) + rotateCalipers(data, graph, values)))
			except ValueError:
				print s + " is not in valid JSON format."

if __name__ == '__main__':
	if len(sys.argv) >= 2:
		args = sys.argv
		main(args[1:])
	else:
		raise LookupError("Please include a schema for which to estimate a cost.")
