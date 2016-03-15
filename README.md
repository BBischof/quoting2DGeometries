# Price quoting 2D geometries of lines and circular arcs

Given a set of straight edges and circular arcs--defined by their center, two points on the arc, and the order of incidence--we generate a minimal rectangular bounding box and calculate the arc length.

Our strategy is similar to the rotating-calipers, i.e. we generate a convex hull of the points, and then iterate over bounding boxes with one edge coincident on one of the hull's edges. For point clouds, the problem is reduced to this method by the theorem of [Freeman and Shapira](http://dl.acm.org/citation.cfm?id=360919), and implemented for this problem in [Arnon and Gieselmann](http://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=1381&context=cstech). 

However, our case differs from the classical. In particular, we need consider extruding curves. We avoid this issue by constructing minimum bounding rectangles around these arcs separately, and unioning these with our convex hull, and then recalculating a new convex hull. It's clear that this produces the minimum convex hull, but it is *not* clear that this will lead to the minimum rectangular bounding box. Some geometries that *may* lead to this problem will be discussed below, and could provide the motivation for further work.

Once the geometric leg-work is done, we use a few constants and a scaling function to calculate a price. 

## Usage

Running the script from the command line takes multiple arguments and returns the price quote for each.

E.g.
```
> python geomCalc.py Rectangle.json
```

Output:
```
Rectangle.json $14.10
```

E.g.
```
> python geomCalc.py Rectangle.json CutCircularArc.json ExtrudeCircularArc.json DoubleRectangle.json lowTop.json
```

Output:
```
Rectangle.json $14.10
CutCircularArc.json $4.06
ExtrudeCircularArc.json $4.47
DoubleRectangle.json $15.48
lowTop.json $14.21
```

### Multiple Closed Curves

Schemas with multiple close curves are supported. They are handled by calculating individual bounding boxes for each(including the padding in each), and combining the cost. This is to reflect that they can be cut individually. 

Note that it might be possible to achieve slightly lower quotes if they were to be cut together. Another approach would be to attempt to redesign the layout such that they are as close to one another as possible before quoting. This could further minimize material usage. In all these cases, arc length is invariant. 

### Special Errors

Failure to provide any argument:
```
"Please include a schema for which to estimate a cost."
```

Failure to provide a valid JSON file:
```
"(arg) is not in valid JSON format.""
```

Failure to provide a JSON in the appropriate Schema format:
```
"(arg) is not a valid schema"
```

Failure to provide a schema consisting of closed curves:
```
"No closed cycle found."
```
Note that if your schema consists of one closed curve and some number of unclosed curves, the closed curve will be quoted, and the others ignored. This functionality could be modified to notify when any extraneous edges are included.

### Outside Sources

The main package utilized is [networkx](https://networkx.github.io/documentation/latest/index.html). The simplicity and guarenteed efficiency was valued over a homebrewed graph implementation. Additionally, the implementation of Graham's convex hull algorithm, [here](https://gist.github.com/tixxit/242402), was used. Again, simplicity and efficiency of the implementation was prioritized.

## Definition Of Schemas

We consider Schemas of the form:
```
{
    "Edges": [
        id: {
            "Type": "LineSegment",
            "Vertices": [id],
        },
        id: {
            "Type": "CircularArc",
            "Center": {
                "X": double,
                "Y": double,
            },
            "ClockwiseFrom": id,
            "Vertices": [id],
        }
    ],
    "Vertices": [
        id: {
            "Position": {
                "X": double,
                "Y": double,
            }
        }
    ]
}
```

## Potentially Dangerous Geometries

It is not clear that smaller calipers may not be found, than those necessary to contain the `extended hulls`. For example, if a circular arc's chord was interior to an `initial hull`'s edge, but the top of the arc extruded beyond this edge, such that the chord and the hull's edge were not parallel, then it is possible that a different, smaller bounding box--than the obvious one calculated using sector height--may be possible to contain the arc. Roughly, one would need to find the point on the arc at the greatest distance from the hull's edge. This point can be found using planar geometry--in particular, for a certain domain of angles between the edge and the chord, this can be obtained by finding the point where the line perpendicular to the edge and intersecting the segment's midpoint, intersects the circular arc. This was not implemented in the final code due to the difficulty in determining the domain of angles and converting this observation to code. A prototypical implementation showed errors and was abandoned; however, given more time, could be rectified.

Note, it's possible that it's actually necessary to calculate these distances from the initial hull edges. Another approach would be to forego computing a secondary hull, and instead involve these arguments in the calculation of the minimum bounding rectangles. This was briefly explored and abandoned for similar reasons to above.

## Functions For Avoiding Graphs

For simplification purposes, I used graphs to represent the connection topology of the schema. This was not necessary, but surely simplified some functions. If one would like to eschew the use of graphs, we can write more tedious implementations of some of these functions. For example:

```
def findPath(name, schema, edges): #attempts to find a path through the edges, keeps for later
	edge_num = len(edges)
	paths = []
	path = [(edges[0], [x  for x in 
		sorted(schema["Edges"][edges[0]]["Vertices"])])]
	tail = path[-1][1][1] #the vertex I'm currently trying to connect up with
	iterations = 0
	while len([s for s in path if tail in s[1]]) == 1: #since I start with the second point in the vertex list, if I see two copies of the first, I've cycled.
		temp = [((v, sorted(schema["Edges"][v]["Vertices"])), #the next edge and its vertices
			(sorted(schema["Edges"][v]["Vertices"]).index(tail)+1)%2) #this is a little trick to find the index of the new tail point
			for v in edges  
			if (tail in schema["Edges"][v]["Vertices"] and v not in [n[0] for n in path])] #grabs another edge that isn't in the path yet, but is incident on the tail
		tail = temp[0][0][1][temp[0][1]]
		path.append(temp[0][0])
		iterations +=1
		if iterations > edge_num:
			raise ValueError(name + " seems to have an unclosed path.")
	paths.append(path)
	if len(path) < edge_num: #recursive call for disconnected paths
		paths += findPath(name, schema, [e for e in edges if e not in [n[0] for n in path]])
	return paths
```

and 

```
def computeInnerHull(name, schema): #helper function to define the first convex hull before extrusions are considered
	hull_points = []
	for path in findPath(name, schema, schema["Edges"].keys()): #compute a hull for each path
		point_names = [] 
		point_names.extend(p for p in itertools.chain.from_iterable([e[1] for e in path]) 
			if p not in point_names)
		hull_points.append(grahms_hull([(schema["Vertices"][str(p)]["Position"]["X"], schema["Vertices"][str(p)]["Position"]["Y"]) 
			for p in point_names]))
	return hull_points
```
Note that these are primitive implementations not intended for production.

## Further Work

In the future, I'd be curious to explore some more things:
- Finding proof or explicit counterexample that extending the convex hull yields the minimum bounding box.
- Improve my schema checking to try to catch more oddities from the outset, such as intersecting edges, or missing data.
- Clean up the rectangles function; both the structure and complexity.
- Implement a function to optimize schemas with multiple curves.