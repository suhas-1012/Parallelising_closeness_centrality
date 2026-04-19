import os, glob, re

files = glob.glob('0[1-7]*.cpp')
for f in files:
    with open(f, 'r') as file:
        content = file.read()
    
    # We want to keep ONLY the line that starts with something like `cout << "Time: "` or similar.
    # What are the time output lines? Usually: cout << "Time: " << ... ms << endl;
    # Let's just comment out all `cout <<` that don't have "Time" or "ms" or "closeness" etc.
    # Actually, it's safer to just comment out the Top 10 block and Method/Formula lines.
    
    # Find Method/Formula/Graph lines:
    content = re.sub(r'(cout\s*<<\s*"Method:.*?;)', r'// \1', content)
    content = re.sub(r'(cout\s*<<\s*"Formula:.*?;)', r'// \1', content)
    content = re.sub(r'(cout\s*<<\s*"Graph:.*?;)', r'// \1', content)
    content = re.sub(r'(cout\s*<<\s*"BCCs:.*?;)', r'// \1', content)
    content = re.sub(r'(cout\s*<<\s*"Articulation points:.*?;)', r'// \1', content)
    content = re.sub(r'(cout\s*<<\s*"R3 nodes:.*?;)', r'// \1', content)
    content = re.sub(r'(cout\s*<<\s*"R4 nodes:.*?;)', r'// \1', content)
    content = re.sub(r'(cout\s*<<\s*"Total redundant:.*?;)', r'// \1', content)
    content = re.sub(r'(cout\s*<<\s*"BFS calls saved:.*?;)', r'// \1', content)
    
    # Top 10 block:
    # cout << "Top 10:" << endl;
    # for (int i = 0; i < min(10, g.n); i++) {
    #     cout << "  Node " << idx[i] << ": " << fixed << setprecision(64) << cc[idx[i]] << endl;
    # }
    content = re.sub(r'(cout\s*<<\s*"Top 10.*?;)', r'// \1', content)
    content = re.sub(r'(cout\s*<<\s*"  Node ".*?;)', r'// \1', content)
    
    with open(f, 'w') as file:
        file.write(content)

