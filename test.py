import subprocess, os, sys
count = 100000
calculate1 = 0.0
calculate2 = 0.0
calculate3 = 0.0
for i in range(0,count):   
    process = subprocess.Popen(['./keygen','sk.key','pk.key'],
                    stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    string = stdout.decode("utf-8")
    sstring = string.split(":")    
    calculate1 += float(sstring[1])


# os.system("echo keygen : "+str(calculate/count)+">out.txt")

# calculate = 0.0
# for i in range(0,count):
    process = subprocess.Popen(['./sign','sk.key','test.txt','sig.txt'],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    string = stdout.decode("utf-8")
    sstring = string.split(":")
    calculate2 += float(sstring[1])



#os.system("echo sign : "+str(calculate/count)+">>out.txt")

# calculate = 0.0
# for i in range(0,count):
    process3 = subprocess.Popen(['./verify','pk.key','test.txt','sig.txt'],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
    stdout, stderr = process3.communicate()
    # os.system("echo "+str(i)+">out.txt")
    string3 = stdout.decode("utf-8")
    sstring3 = string3.split(":")
    calculate3 += float(sstring3[1])
print("Keygen:",calculate1/count)
print("Sign:",calculate2/count)
print("Verify:",calculate3/count)

# os.system("echo verify : "+str(calculate/count)+">out.txt")
