import re

files = ['08_novel_bcc_spmm.cpp', '09_novel_vdbcc.cpp']

for f in files:
    with open(f, 'r') as file:
        content = file.read()
    
    # Remove bfs_local completely
    content = re.sub(r'//bfs within a local bcc\nstatic vector<int> bfs_local.*?\n}\n\n', '\n', content, flags=re.DOTALL)
    content = re.sub(r'/\* BFS within a local BCC \*/\nstatic vector<int> bfs_local.*?\n}\n\n', '\n', content, flags=re.DOTALL)
    
    # Remove dynamicInsert support
    content = re.sub(r'//dynamic edge insertion support\nstruct DynamicResult.*?\n}\n\n', '\n', content, flags=re.DOTALL)
    content = re.sub(r'/\* ────── Dynamic Edge Insertion Support ────── \*/\nstruct DynamicResult.*?\n}\n\n', '\n', content, flags=re.DOTALL)
    
    # Remove baseline_cc
    content = re.sub(r'//baseline sequential BFS CC\nvoid baseline_cc.*?\n}\n\n\n', '\n', content, flags=re.DOTALL)
    content = re.sub(r'/\* ────── Baseline sequential BFS CC ────── \*/\nvoid baseline_cc.*?\n}\n\n', '\n', content, flags=re.DOTALL)
    
    with open(f, 'w') as file:
        file.write(content)
