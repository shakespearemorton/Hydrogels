for i in range(25):
    str1 = str('system' +repr(i))
    with open('template.in', 'r') as file :
        filedata = file.read()
    filedata = filedata.replace('REPLACE', str1)
    with open('template_'+repr(i)+'.in', 'w') as file:
        file.write(filedata)
        
    with open('run.pbs', 'r') as file :
        filedata = file.read()

    filedata = filedata.replace('REPLACE', repr(i))
    
    with open('run_'+repr(i)+'.pbs', 'w') as file:
        file.write(filedata)