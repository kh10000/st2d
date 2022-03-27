import pickle, sys

fil = sys.argv[1]

with open(fil, "rb") as f:
    res = pickle.load(f)

print("Sr")
print(res["x"]["Sr"])
print("Ti")
print(res["x"]["Ti"])
print("O")
print(res["x"]["O"])