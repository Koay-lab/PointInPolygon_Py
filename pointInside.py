# the representation of a point will be a tuple (x,y)
# the representation of a polygon wil be a list of points [(x1,y1), (x2,y2), (x3,y3), ... ]

# it is assumed that polygon is regular i.e. lines don't intersect each other (otherwise, it is questionable whether it is a polygon)

def isLeft(P0, P1, P2):
    """
    Tests if a point is Left|On|Right of an infinite line.
    Input:  three points P0, P1, and P2
    :return: >0 for P2 left of the line through P0 and P1, =0 for P2  on the line, <0 for P2  right of the line
    """
    return (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P2[0] - P0[0]) * (P1[1] - P0[1])


def wn_PnPoly(P, V):
    """
    Winding number test for a point in a polygon.
    Taken from: http://geomalgorithms.com/a03-_inclusion.html
    :param P: a point
    :param V: vertex points of a polygon V[n+1] with V[n]=V[0]
    :return: wn = the winding number (=0 only when P is outside)
    """
    n = len(V) - 1
    wn = 0  # the  winding number counter
    i = 0

    # loop through all edges of the polygon
    while i < n:  # edge from V[i] to  V[i+1]
        if V[i][1] <= P[1]:  # start y <= P.y
            if V[i + 1][1] > P[1]:  # an upward crossing
                if isLeft(V[i], V[i + 1], P) > 0:  # P left of  edge
                    wn += 1  # have  a valid up intersect
        else:  # start y > P.y (no test needed)
            if V[i + 1][1] <= P[1]:  # a downward crossing
                if isLeft(V[i], V[i + 1], P) < 0:  # P right of  edge
                    wn -= 1  # have  a valid down intersect

        i += 1

    # print str(wn)
    return wn


def cn_PnPoly(P, V):  # P - Point; V - Polygon
    """
    Crossing number test for a point in a polygon.
    Taken from: http://geomalgorithms.com/a03-_inclusion.html
    This code is patterned after [Franklin, 2000]

    :param P: a point,
    :param V: vertex points of a polygon V[n+1] with V[n]=V[0]
    :return: 0 = outside, 1 = inside
    """
    n = len(V) - 1
    cn = 0  # the  crossing number counter

    # loop through all edges of the polygon
    i = 0
    # while i<n: # edge from V[i]  to V[i+1]
    for i in range(n):
        # upward crossing or downward crossing
        if (V[i][1] <= P[1] and V[i + 1][1] > P[1]) or (V[i][1] > P[1] and V[i + 1][1] <= P[1]):
            # compute  the actual edge-ray intersect x-coordinate
            vt = (P[1] - V[i][1]) / (V[i + 1][1] - V[i][1])
            # if P[0] > (V[i][0] + vt * (V[i+1][0] - V[i][0])): # P.x > intersect - ray toward left
            if P[0] < (V[i][0] + vt * (V[i + 1][0] - V[i][0])):  # P.x < intersect - ray toward right - original
                cn += 1  # a valid crossing of y=P[1] right of P.x

    # print str(cn)
    return (cn & 1)  # 0 if even (out), and 1 if  odd (in)


# -------------------------------------------------------------------------------
def isPointInPath(P, poly):
    """
    Taken from: https://en.wikipedia.org/wiki/Even%E2%80%93odd_rule

    :param P: a point
    :param poly: vertex points of a polygon poly[n+1] with poly[n]=poly[0]
    :return: False = outside, True = inside
    """
    num = len(poly)
    i = 0
    j = num - 1
    c = False
    for i in range(num):
        if ((poly[i][1] > P[1]) != (poly[j][1] > P[1])) and \
                (P[0] < (poly[j][0] - poly[i][0]) * (P[1] - poly[i][1]) / (poly[j][1] - poly[i][1]) + poly[i][0]):
            c = not c
        j = i

    return c


# -------------------------------------------------------------------------------
def is_inside_postgis(polygon, point):
    length = len(polygon)
    intersections = 0

    dx2 = point[0] - polygon[0][0]
    dy2 = point[1] - polygon[0][1]
    ii = 0
    jj = 1

    while jj < length:
        dx = dx2
        dy = dy2
        dx2 = point[0] - polygon[jj][0]
        dy2 = point[1] - polygon[jj][1]

        F = (dx - dx2) * dy - dx * (dy - dy2);
        if 0.0 == F and dx * dx2 <= 0 and dy * dy2 <= 0:
            return 2;

        if (dy >= 0 and dy2 < 0) or (dy2 >= 0 and dy < 0):
            if F > 0:
                intersections += 1
            elif F < 0:
                intersections -= 1

        ii = jj
        jj += 1

    # print 'intersections =', intersections
    return intersections != 0


# -------------------------------------------------------------------------------
def is_inside_sm(polygon, point, tol=0):
    """
    This function gives the answer whether the given point is inside or outside the predefined polygon.
    Unlike standard ray-casting algorithm, this one works on edges! (with no performance cost)
    According to performance tests - this is the best variant.

    The algorithm is based on counting the number of edges of the polygon that are to the left of and not completely
    above or below the point of interest. If and only if this number is odd, then the point is inside the polygon
    (since it is a closed region).

    :param polygon: searched polygon, as a list of vertex points with polygon[-1] == polygon[0]
    :param point: an arbitrary point that can be inside or outside the polygon
    :return: 0 - the point is outside the polygon; 1 - the point is inside the polygon; 2 - the point is one edge (boundary)
    """

    length = len(polygon) - 1
    dy2 = point[1] - polygon[0][1]
    intersections = 0
    ii = 0
    jj = 1

    # point on endpoint of line (dy2=dx2=0)
    if abs(dy2) <= tol and abs(point[0] - polygon[0][0]) <= tol:
        return 2

    while ii < length:
        dy = dy2
        dy2 = point[1] - polygon[jj][1]

        # plt.plot([polygon[ii][0], polygon[jj][0]], [polygon[ii][1], polygon[jj][1]]); plt.scatter(*point); plt.show()

        # point on endpoint of line (dy2=dx2=0)
        if abs(dy2) <= tol and abs(point[0] - polygon[jj][0]) <= tol:
            return 2

        # horizontal line (dy == dy2)
        if abs(dy - dy2) <= tol:
            # point on horizontal line (dy=dy2=0 and dx*dx2<=0)
            if abs(dy) <= tol and abs(dy2) <= tol and (point[0] - polygon[ii][0]) * (point[0] - polygon[jj][0]) <= tol:
                return 2

        # consider only lines which are not completely above/below/right of the point
        elif ((dy < 0) ^ (dy2 < 0)) and (point[0] + tol >= polygon[ii][0] or point[0] + tol >= polygon[jj][0]):
            F = dy * (polygon[jj][0] - polygon[ii][0]) / (dy - dy2) + polygon[ii][0]

            if abs(point[0] - F) <= tol:  # point on line
                return 2
            if point[0] > F:  # if line is left from the point - the ray moving towards left, will intersect it
                intersections += 1

        # there is another possibility: (dy=0 and dy2>0) or (dy>0 and dy2=0). It is skipped
        # deliberately to prevent break-points intersections to be counted twice.

        ii = jj
        jj += 1

    # print 'intersections =', intersections
    return intersections & 1
