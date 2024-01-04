from quiver import *

def random_test_data(vertices, arrowBound, dimBound):
    inFundamentalDomain = False
    while not inFundamentalDomain:
        Q = RandomQuiver(vertices, arrow_bound=arrowBound, acyclic=True, connected=True)
        d = RandomDimensionVector(Q, positive=True, upper_bound=dimBound)
        inFundamentalDomain = Q.in_fundamental_domain(d)
    return Q, d

vertices = 3
arrowBound = 10
dimBound = 10

sampleSize = 1000
stronglyAmplyStableCount = 0

for i in range(sampleSize):
    Q, d = random_test_data(vertices, arrowBound, dimBound)
    theta = Q.canonical_stability_parameter(d)
    if Q.is_strongly_amply_stable(d, theta):
        stronglyAmplyStableCount += 1
    
print("Test random quivers with "+str(vertices)+" vertices and random dimension vectors for strong ample stability.")
print("Sample size: "+str(sampleSize))
print("Strongly amply stable count: "+str(stronglyAmplyStableCount))
