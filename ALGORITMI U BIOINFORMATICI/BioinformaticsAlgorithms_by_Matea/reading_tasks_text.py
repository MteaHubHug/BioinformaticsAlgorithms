import os

##### 1 A ######
def find_file(name):
    path=r"C:\Users\Matea\Desktop\ALGORITMI U BIOINFORMATICI\BioinformaticsAlgorithms_by_Matea\tekstovi_zadataka"
    text_files=os.listdir(path)
    for file in text_files:
        if(name==file):
            return path+"\\"+file

def open_file(name):
    file=find_file(name)
    f = open(file, "r")
    content=f.readlines()
    content_edited=[]
    for line in content:
        l=line.strip('\n')
        content_edited.append(l)
    return content_edited


