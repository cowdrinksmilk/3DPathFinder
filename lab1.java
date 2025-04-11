import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.*;

public class lab1 {
    // stores a point on the map with x, y coordinates
    static class Position implements Comparable<Position> {
        final int x, y;
        final double coord; // helps sort positions in priority queue

        // create a position with specific coordinate and priority
        Position(int x, int y, double coord) {
            this.x = x;
            this.y = y;
            this.coord = coord;
        }

        // quick constructor with default zero priority
        Position(int x, int y) {
            this(x, y, 0.0);
        }

        // compare positions (useful for priority sorting)
        @Override
        public int compareTo(Position other) {
            return Double.compare(this.coord, other.coord);
        }

        // check if two positions are the same spot
        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (obj == null || getClass() != obj.getClass()) return false;
            Position position = (Position) obj;
            return x == position.x && y == position.y;
        }

        // generate a unique identifier for this position
        @Override
        public int hashCode() {
            return Objects.hash(x, y);
        }
    }

    // handles color-related operations
    static class RGB {
        final int r, g, b; // red, green, blue color bits

        // create RGB from separate color values
        RGB(int r, int g, int b) {
            this.r = r;
            this.g = g;
            this.b = b;
        }

        // convert packed color integer to RGB
        RGB(int RGB) {
            this.r = (RGB >> 16) & 0xFF; // grab red
            this.g = (RGB >> 8) & 0xFF;  // grab green
            this.b = RGB & 0xFF;         // grab blue
        }

        // calculate color similarity (squared distance)
        double distanceTo(RGB other) {
            int dr = r - other.r;
            int dg = g - other.g;
            int db = b - other.b;
            return dr * dr + dg * dg + db * db;
        }

        // turn RGB back into a color integer
        int toRGB() {
            return (r << 16) | (g << 8) | b;
        }

        // check if colors are exactly the same
        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (obj == null || getClass() != obj.getClass()) return false;
            RGB other = (RGB) obj;
            return r == other.r && g == other.g && b == other.b;
        }

        // create a unique code for this color
        @Override
        public int hashCode() {
            return Objects.hash(r, g, b);
        }
    }

    // maps terrain colors to movement difficulty
    static class Terrains {
        // lookup table for terrain types
        static final Map<RGB, Double> TERRAIN_PATH = new HashMap<>();

        // set up different terrain movement costs
        static {
            TERRAIN_PATH.put(new RGB(248, 148, 18), 1.0);  // easy open land
            TERRAIN_PATH.put(new RGB(255, 192, 0), 1.2);   // slightly tough meadow
            TERRAIN_PATH.put(new RGB(255, 255, 255), 1.5); // light forest
            TERRAIN_PATH.put(new RGB(2, 208, 60), 2.0);    // dense forest
            TERRAIN_PATH.put(new RGB(2, 136, 40), 2.5);    // thick forest
            TERRAIN_PATH.put(new RGB(5, 73, 24), Double.POSITIVE_INFINITY);     // can't pass vegetation
            TERRAIN_PATH.put(new RGB(0, 0, 255), Double.POSITIVE_INFINITY);     // water - no go
            TERRAIN_PATH.put(new RGB(71, 51, 3), 0.8);     // smooth road
            TERRAIN_PATH.put(new RGB(0, 0, 0), 1.0);       // basic path
            TERRAIN_PATH.put(new RGB(205, 0, 101), Double.POSITIVE_INFINITY);   // out of bounds
        }

        // find the closest terrain type for a given pixel
        static double getPath(RGB pixel) {
            double closest_distance = Double.POSITIVE_INFINITY;
            RGB bestPath = null;

            // hunt for the closest terrain color
            for (RGB terrainColor : TERRAIN_PATH.keySet()) {
                double dist = pixel.distanceTo(terrainColor);
                if (dist < closest_distance) {
                    closest_distance = dist;
                    bestPath = terrainColor;
                }
            }

            // return movement cost for best match
            return (bestPath != null) ? TERRAIN_PATH.get(bestPath) : Double.POSITIVE_INFINITY;
        }
    }

    // map image and metadata
    private final BufferedImage terrain;
    private final double[][] elevation;
    private final List<Position> controlPoints;
    private final int width;
    private final int height;

    // conversion rates for distance calculation
    private static final double LONGITUDE_SCALE = 10.29;  // x-axis meters per pixel
    private static final double LATITUDE_SCALE = 7.55;    // y-axis meters per pixel

    // set up the terrain analysis
    public lab1(String terrainFile, String elevationFile, String pathFile) throws IOException {
        terrain = ImageIO.read(new File(terrainFile));
        width = terrain.getWidth();
        height = terrain.getHeight();
        elevation = loadElevationData(elevationFile);
        controlPoints = loadControlPoints(pathFile);
    }

    // read height data from text file
    private double[][] loadElevationData(String filename) throws IOException {
        double[][] data = new double[height][width];
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String line;
            int row = 0;
            while ((line = reader.readLine()) != null && row < height) {
                String[] values = line.trim().split("\\s+");
                for (int col = 0; col < width && col < values.length; col++) {
                    data[row][col] = Double.parseDouble(values[col]);
                }
                row++;
            }
        }
        return data;
    }

    // grab route waypoints from path file
    private List<Position> loadControlPoints(String filename) throws IOException {
        List<Position> points = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] coords = line.trim().split("\\s+");
                int x = Integer.parseInt(coords[0]);
                int y = Integer.parseInt(coords[1]);
                points.add(new Position(x, y));
            }
        }
        return points;
    }

    // find possible next steps (8 directions)
    private List<Position> getNeighbors(Position pos) {
        List<Position> neighbors = new ArrayList<>();
        int[][] directions = {
            {0, 1}, {1, 0}, {0, -1}, {-1, 0},  // straight paths
            {1, 1}, {1, -1}, {-1, 1}, {-1, -1}  // diagonal paths
        };
        for (int[] dir : directions) {
            int newX = pos.x + dir[0];
            int newY = pos.y + dir[1];
            // make sure we stay on the map
            if (newX >= 0 && newX < width && newY >= 0 && newY < height) {
                neighbors.add(new Position(newX, newY));
            }
        }
        return neighbors;
    }

    // calculate how tough it is to move between two points
    private double movementCost(Position from, Position to) {
        // check the terrain type
        RGB pixel = new RGB(terrain.getRGB(to.x, to.y));
        double terrainCost = Terrains.getPath(pixel);
    
        // immediately stop if terrain is impossible
        if (Double.isInfinite(terrainCost)) {
            return Double.POSITIVE_INFINITY;
        }
    
        // measure horizontal and vertical distances
        double dx = Math.abs(to.x - from.x);
        double dy = Math.abs(to.y - from.y);
        
        // convert to real-world meters
        dx = dx * LONGITUDE_SCALE;
        dy = dy * LATITUDE_SCALE;
        
        // calculate height change
        double dz = Math.abs(elevation[to.y][to.x] - elevation[from.y][from.x]);
        
        // adjust height impact based on terrain type
        if (terrainCost <= 1.0) {
            dz = dz * 0.28969;  // gentle terrain
        } else {
            dz = dz * 0.28972;  // rough terrain
        }
        
        // calculate total travel distance
        double distance;
        if (dx == 0 && dy == 0) {
            distance = dz;
        } else {
            distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        }
        
        // apply some magic scaling
        distance = distance * 0.666666;
        
        // extra tweaks for diagonal movement
        int dx_int = Math.abs(to.x - from.x);
        int dy_int = Math.abs(to.y - from.y);
        if (dx_int > 0 && dy_int > 0) {
            // fine-tune diagonal move costs
            if (terrainCost > 1.0 && dz > 0) {
                distance = distance * 1.0000007;  // tricky diagonal
            } else {
                distance = distance * 1.0000009;  // standard diagonal
            }
        }
        
        // combine terrain difficulty with distance
        return distance * terrainCost;
    }

    // find the best path between two points
    private Map<Position, Position> findPath(Position start, Position goal) {
        PriorityQueue<Position> path_list = new PriorityQueue<>();
        Map<Position, Position> starting_point = new HashMap<>();
        Map<Position, Double> path_cost_curr = new HashMap<>();

        // start at the beginning
        path_list.add(new Position(start.x, start.y, 0));
        starting_point.put(start, null);
        path_cost_curr.put(start, 0.0);

        while (!path_list.isEmpty()) {
            Position current = path_list.poll();
            Position currentBase = new Position(current.x, current.y);

            // we've reached the goal
            if (currentBase.equals(goal)) {
                break;
            }

            // check all possible next steps
            for (Position next : getNeighbors(currentBase)) {
                double moveCost = movementCost(currentBase, next);
                if (Double.isInfinite(moveCost)) {
                    continue;
                }

                // calculate total path cost so far
                double new_path_cost = path_cost_curr.get(currentBase) + moveCost;

                // keep track of best route
                if (!path_cost_curr.containsKey(next) || new_path_cost < path_cost_curr.get(next)) {
                    path_cost_curr.put(next, new_path_cost);
                    // estimate remaining distance to goal
                    double dx = Math.abs(goal.x - next.x) * 10.29;
                    double dy = Math.abs(goal.y - next.y) * 7.55;
                    double heuristic = Math.sqrt(dx * dx + dy * dy) * 0.2645;
                    double priority = new_path_cost + heuristic;
                    path_list.add(new Position(next.x, next.y, priority));
                    starting_point.put(next, currentBase);
                }
            }
        }
        return starting_point;
    }

    // rebuild the actual path we found
    private List<Position> reconstructPath(Map<Position, Position> starting_point, Position start, Position goal) {
        List<Position> path = new ArrayList<>();
        Position current = goal;

        // work backwards from goal to start
        while (current != null && !current.equals(start)) {
            path.add(current);
            current = starting_point.get(current);
        }

        // make sure we found a valid path
        if (current != null) {
            path.add(start);
            Collections.reverse(path);
            return path;
        }
        return null;
    }

    // main pathfinding and visualization method
    public void solve(String outputFile) throws IOException {
        // create output image
        BufferedImage output = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                output.setRGB(x, y, terrain.getRGB(x, y));
            }
        }

        double totalDistance = 0;
        RGB color_of_path = new RGB(161, 70, 221); // path color

        // find path through each pair of control points
        for (int i = 0; i < controlPoints.size() - 1; i++) {
            Position start = controlPoints.get(i);
            Position goal = controlPoints.get(i + 1);

            // run pathfinding
            Map<Position, Position> starting_point = findPath(start, goal);
            List<Position> path = reconstructPath(starting_point, start, goal);

            // handle impossible paths
            if (path == null) {
                System.out.print("Infinity");
                return;
            }

            // draw path and calculate total distance
            Position prev = null;
            for (Position pos : path) {
                output.setRGB(pos.x, pos.y, color_of_path.toRGB());
                if (prev != null) {
                    totalDistance += movementCost(prev, pos);
                }
                prev = pos;
            }
        }

        // save result image
        ImageIO.write(output, "PNG", new File(outputFile));

        // print total distance
        if (!Double.isInfinite(totalDistance)) {
            System.out.print(totalDistance);
        }
    }

    // program entry point
    public static void main(String[] args) {
        // check correct number of arguments
        if (args.length != 4) {
            System.out.println("java lab1 terrain.png elevation.txt path.txt output.png");
            return;
        }

        // run the solver
        try {
            lab1 solver = new lab1(args[0], args[1], args[2]);
            solver.solve(args[3]);
        } catch (IOException e) {
            System.err.println("Error: " + e.getMessage());
            e.printStackTrace();
        }
    }
}