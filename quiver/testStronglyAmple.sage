from quiver import *

def random_test_data(vertices, arrowBound, dimBound):
    inFundamentalDomain = False
    while not inFundamentalDomain:
        Q = RandomQuiver(vertices, arrow_bound=arrowBound, acyclic=True, connected=True)
        d = RandomDimensionVector(Q, positive=True, upper_bound=dimBound)
        inFundamentalDomain = Q.in_fundamental_domain(d)
    return Q, d

print("Test random Q and d for strong ample stability, where d is in fundamental domain.")
print("Number of vertices: ")
vertices = int(input())
arrowBound = 10
dimBound = 10

print("Sample size: ")
sampleSize = int(input())
stronglyAmplyStableCount = 0

for i in range(sampleSize):
    Q, d = random_test_data(vertices, arrowBound, dimBound)
    theta = Q.canonical_stability_parameter(d)
    if Q.is_strongly_amply_stable(d, theta):
        stronglyAmplyStableCount += 1
    
print("\nNumber of vertices: "+str(vertices))
print("Upper bound for number of arrows between two vertices: "+str(arrowBound))
print("Upper bound for entries of d: "+str(dimBound)+"\n")
print("Sample size: "+str(sampleSize))
print("Strongly amply stable count: "+str(stronglyAmplyStableCount))
