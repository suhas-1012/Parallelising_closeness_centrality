import glob, re

for f in glob.glob('0[1-7]*.cpp'):
    with open(f, 'r') as file:
        content = file.read()
    
    # Remove the `for (int i = 0; i < min(10...` entirely or comment it out
    content = re.sub(r'(for\s*\([^\{]*min\(10.*?g\.n\)[^\{]*\))', r'// \1', content)
    
    with open(f, 'w') as file:
        file.write(content)
