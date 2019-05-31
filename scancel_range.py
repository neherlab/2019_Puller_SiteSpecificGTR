import sys, os


for i in range(int(sys.argv[1]), int(sys.argv[2])):
    os.system("scancel " +str(i))
