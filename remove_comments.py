import re

def remove_comments(text):
    pattern = r"//.*?$|/\*.*?\*/|'(?:\\.|[^\\'])*'|\"(?:\\.|[^\\\"])*\""
    def replacer(match):
        s = match.group(0)
        if s.startswith('/'):
            return "" 
        else:
            return s
    return re.sub(pattern, replacer, text, flags=re.DOTALL | re.MULTILINE)

with open("09_novel_vdbcc.cpp", "r") as f:
    content = f.read()

content = remove_comments(content)
# collapse more than 2 blank lines into 2
content = re.sub(r'\n{3,}', '\n\n', content)

with open("09_novel_vdbcc.cpp", "w") as f:
    f.write(content)
