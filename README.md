# Orienteering Path Finder

This project solves a digital orienteering problem using A* search. Given a terrain map, elevation data, and a list of control points to visit in order, it finds the fastest path across the map while accounting for both terrain types and elevation changes.

## What it does

- Calculates the shortest-time path between a sequence of points  
- Uses color-coded terrain and elevation to determine movement cost  
- Draws the final path on the original map image  
- Prints the total path length (in meters) to the terminal

## How to run it

Make sure you have Python 3 installed, and install dependencies:

```bash
pip install pillow numpy
```

Then run the program like this:

```bash
python3 lab1.py terrain.png elevation.txt controls.txt output.png
```

- `terrain.png`: the map image (395x500, color-coded by terrain type)  
- `elevation.txt`: elevation values (500 rows of 400 floats — ignore last 5 values per row)  
- `controls.txt`: list of (x, y) points to visit  
- `output.png`: name for the output image with the path drawn on it

## Output

- A new image file with the path drawn in purple (`#a146dd`)  
- The total distance of the path printed to the terminal

## Notes

- Movement is affected by terrain type (e.g., forest is slower than roads)
- Elevation is used in the 3D distance calculation
- The control points must be visited in order
- The path is drawn 1 pixel wide

---
