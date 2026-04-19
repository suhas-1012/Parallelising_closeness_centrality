import glob
for f in glob.glob('0[1-7]*.cpp'):
    with open(f, 'r') as file:
        content = file.read()
    
    content = content.replace('for (int i = 0; i < min(10, g.n); i++)', '// for (int i = 0; i < min(10, g.n); i++)')
    content = content.replace('for (int i = 0; i < min(10, n); i++)', '// for (int i = 0; i < min(10, n); i++)')
    
    with open(f, 'w') as file:
        file.write(content)
