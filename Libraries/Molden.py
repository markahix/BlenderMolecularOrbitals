def extract_coordinate_lines(molden):
    with open(molden,"r") as f:
        lines = f.readlines()
    startprinting = False
    coord_lines = []
    for line in lines:
        if "[GTO]" in line:
            startprinting = False
        if startprinting:
            coord_lines.append(line)
        if "[Atoms] Angs" in line:
            startprinting = True

    return coord_lines
