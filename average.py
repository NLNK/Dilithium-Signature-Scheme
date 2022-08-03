import subprocess, os, sys
count = 100000
calculate = 0.0
for i in range(0,count):   
    process = subprocess.Popen(['./keygen','sk.key','pk.key'],
                    stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    string = stdout.decode("utf-8")
    sstring = string.split(":")    
    calculate += float(sstring[1])

print("Keygen:",calculate/count)
# os.system("echo keygen : "+str(calculate/count)+">out.txt")

calculate = 0.0
for i in range(0,count):
    process = subprocess.Popen(['./sign','sk.key','test.txt','sig.txt'],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    string = stdout.decode("utf-8")
    sstring = string.split(":")
    calculate += float(sstring[1])

print("Sign:",calculate/count)

#os.system("echo sign : "+str(calculate/count)+">>out.txt")

calculate = 0.0
for i in range(0,count):
    process = subprocess.Popen(['./verify','pk.key','test.txt','sig.txt'],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    # os.system("echo "+str(i)+">out.txt")
    string = stdout.decode("utf-8")
    sstring = string.split(":")
    calculate += float(sstring[1])

print("Verify:",calculate/count)

# os.system("echo verify : "+str(calculate/count)+">out.txt")
