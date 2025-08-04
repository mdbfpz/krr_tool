from collections import OrderedDict
import asyncio
from geographiclib.geodesic import Geodesic
from shapely.geometry import Point, Polygon, LineString
import math
from datetime import datetime
import matplotlib.pyplot as plt
from shapely.ops import transform
import pyproj
import contextlib
import io

class DataQueue:
    """Manages the queue for storing data asynchronously."""
    
    def __init__(self):
        self.queue = asyncio.Queue()
        self.condition = asyncio.Condition()

    async def add_to_queue(self, data: dict): # TODO: rename to 'put'
        """Add data to the queue asynchronously."""

        async with self.condition:
            await self.queue.put(data)
            self.condition.notify()  # Notify any waiting tasks

    async def get_from_queue(self) -> dict:  # TODO: rename to 'get'
        """Retrieve data from the queue asynchronously."""

        async with self.condition:
            while self.queue.empty():
                await self.condition.wait()  # Wait for new data
            return await self.queue.get()


class GeodesicService:
    """Service class for geodesic calculations."""

    a = 6378137.0  # semi-major axis
    b = 6356752.314245  # semi-minor axis
    c = 0.514444  # Conversion factor for knots to m/s

    @staticmethod
    def check_if_point_in_circle(lat1, lon1, lat2, lon2, radius):
        """
        Check if a point is within a circle defined by its center and radius.
        """
        
        distance = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)['s12']
        return distance <= radius

    @staticmethod
    def calculate_tolerance_azi(lat1, lon1, lat2, lon2, radius):
        """
        Calculate the azimuth tolerance based on radius and distance to the circle center.
        """

        distance = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)['s12']
        if distance != 0:
            alfa_in_radians = math.atan(radius/distance)
            return math.degrees(alfa_in_radians)
        else:
            return 0
    @staticmethod
    def xyz_from_lat_lon(lat, lon, height):
        """
        Convert latitude, longitude, and height to Cartesian coordinates.
        Expects lat and lon in decimal degrees.
        """

        # Convert lat and lon from degrees to radians
        lat_rad = math.radians(lat)
        lon_rad = math.radians(lon)
        print(f"Lat (rad): {lat_rad}, Lon (rad): {lon_rad}")

        cos_lat = math.cos(lat_rad)
        sin_lat = math.sin(lat_rad)
        print(f"Cos(lat): {cos_lat}, Sin(lat): {sin_lat}")

        # Calculate radius of curvature in the prime vertical
        r_n = (GeodesicService.a**2) / math.sqrt(
            (GeodesicService.a**2) * (cos_lat**2) + (GeodesicService.b**2) * (sin_lat**2)
        )
        print(f"Radius of curvature (r_n): {r_n}")

        # Calculate Cartesian coordinates
        x = (r_n + height) * cos_lat * math.cos(lon_rad)
        y = (r_n + height) * cos_lat * math.sin(lon_rad)
        z = ((GeodesicService.b**2 / GeodesicService.a**2) * r_n + height) * sin_lat
        print(f"Computed Cartesian coordinates: x={x}, y={y}, z={z}")

        return x, y, z

    @staticmethod
    def xyz_to_lat_lon(x, y, z):
        """
        Convert Cartesian coordinates to decimal degrees latitude, longitude, and height.
        """

        # Compute the eccentricity squared
        e_2 = (GeodesicService.a**2 - GeodesicService.b**2) / GeodesicService.a**2
        ep_2 = (GeodesicService.a**2 - GeodesicService.b**2) / GeodesicService.b**2
        print(f"e_2: {e_2}, ep_2: {ep_2}")

        # Compute the distance from the origin in the x-y plane
        p = math.sqrt(x**2 + y**2)
        print(f"p (distance in x-y plane): {p}")

        # Compute the angle theta
        theta = math.atan2(z * GeodesicService.a, p * GeodesicService.b)
        print(f"theta: {theta}")

        # Latitude and Longitude in radians
        lon_rad = math.atan2(y, x)
        lat_rad = math.atan2(
            z + ep_2 * GeodesicService.b * math.sin(theta)**3,
            p - e_2 * GeodesicService.a * math.cos(theta)**3
        )
        print(f"Longitude (rad): {lon_rad}, Latitude (rad): {lat_rad}")

        # Convert radians to decimal degrees
        lon = math.degrees(lon_rad)
        lat = math.degrees(lat_rad)
        print(f"Longitude (deg): {lon}, Latitude (deg): {lat}")

        # Height
        r_n = GeodesicService.a / math.sqrt(1 - e_2 * math.sin(lat_rad)**2)
        print(f"Radius of curvature (r_n): {r_n}")
        height = p / math.cos(lat_rad) - r_n
        print(f"Height: {height}")

        return lat, lon, height


    @staticmethod
    def check_if_point_in_polygon(lat, lon, coords):
        """
        Check if a point is inside a polygon (without the edge). 
        Automatically converts polygon coordinates from (lat, lon) to (lon, lat) if necessary.
        """

        polygon, _ = GeodesicService.create_polygon_from_coords(coords)
        point = Point(lon, lat)  # (lon, lat) format

        # Plot
        """print("Point: ", Point(point))
        print("Polygon: ", polygon, "Is valid: ", polygon.is_valid)
        print("Polygon contains point: ", polygon.contains(Point(point)))
        print("Polygon contains point properly: ", polygon.contains_properly(Point(point)))
        print("Exterior: ", polygon.exterior)
        print("Exterior contains point: ", polygon.exterior.contains(Point(point)))
        print("Exterior contains point properly: ", polygon.exterior.contains_properly(Point(point)))
        print("\n")
        GeodesicService.plot_sector_with_annotations(polygon, point)"""

        return polygon.contains(point)
    
    @staticmethod
    def plot_sector_with_annotations(polygon, point, midpoints=None, projections=None):
        """
        Plot the polygon, point, midpoints, and projection points with annotations.
        """

        # Plot the polygon
        lon, lat = zip(*polygon.exterior.coords)
        plt.plot(lon, lat, label="Sector")

        # Plot the main point
        plt.scatter(point.x, point.y, color='red', label="Current position")

        # Annotate midpoints
        if midpoints:
            for idx, (mid_lat, mid_lon) in enumerate(midpoints):
                plt.text(mid_lon, mid_lat, str(idx), color='black', fontsize=12, ha='center', va='center')

        # Plot and annotate projections
        if projections:
            for idx, (proj_lat, proj_lon) in enumerate(projections):
                if proj_lat is not None and proj_lon is not None:
                    plt.scatter(proj_lon, proj_lat, color='green', label=f"Projection Point" if idx == 0 else "")

        # Ensure proper aspect ratio
        plt.gca().set_aspect('equal', adjustable='box')
        plt.legend()
        plt.show()

    @staticmethod
    def geodesic_distance(lat1, lon1, lat2, lon2):
        """Calculates geodesic distance in meters."""

        return Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)['s12']

    @staticmethod
    def calculate_distance_to_sector(lat, lon, coords):
        """
        Calculate the minimum geodesic distance from a point to a polygon in meters.
        """

        # TODO: can this be implemented as: create a line from the point to the projection, then access the distance property?
        # TODO: could we use .project(other, normalized=False) from shapely docs?

        polygon, sorted_cords = GeodesicService.create_polygon_from_coords(coords)
        point = Point(lon, lat)  # (lon, lat) format
        # Plot
        # GeodesicService.plot_sector_with_annotations(polygon, point)
        
        min_distance = float('inf')
        midpoints, projections = [], []
        # Iterate through each edge of the polygon
        for i in range(len(sorted_cords)):
            start = sorted_cords[i]
            end = sorted_cords[(i + 1) % len(sorted_cords)]  # Wrap around to the first vertex

            # Handle degenerate case where start and end are the same
            if start == end:
                # Calculate the distance to the single point (start == end)
                start_distance = GeodesicService.geodesic_distance(lat, lon, start[0], start[1])
                min_distance = min(min_distance, start_distance)
                continue  # Skip to the next edge

            # Calculate the distance to the start and end points (vertices)
            start_distance = GeodesicService.geodesic_distance(lat, lon, start[0], start[1])
            end_distance = GeodesicService.geodesic_distance(lat, lon, end[0], end[1])

            # Calculate the vector from start to end of the segment
            dx = end[1] - start[1]  # Longitude difference
            dy = end[0] - start[0]  # Latitude difference

            # If dx and dy are both zero, skip the projection and continue with the next edge
            if dx == 0 and dy == 0:
                continue  # Degenerate case where the segment is a point

            # Calculate the vector from the start point to the given point (lat, lon)
            px = lon - start[1]
            py = lat - start[0]

            # Calculate the dot product of the vectors (start->end) and (start->point)
            dot_product = px*dx + py*dy

            # Calculate the squared length of the segment (start->end)
            segment_length_squared = dx**2 + dy**2

            # Calculate the projection scalar t
            t = dot_product/segment_length_squared

            if t < 0:
                # The projection falls before the segment, use the start point
                proj_lat, proj_lon = start
            elif t > 1:
                # The projection falls after the segment, use the end point
                proj_lat, proj_lon = end
            else:
                # The projection falls within the segment
                proj_lat = start[0] + t*dy
                proj_lon = start[1] + t*dx
            
            projections.append((proj_lat, proj_lon))
            
            # Calculate the geodesic distance from the point to the projected point
            proj_distance = GeodesicService.geodesic_distance(lat, lon, proj_lat, proj_lon)
            # print(i, ": proj_distance:", proj_distance)

            # Update the minimum distance
            min_distance = min(min_distance, start_distance, end_distance, proj_distance)
            
            # Add the index of the segment (i) at the midpoint of the segment
            mid_lat = (start[0] + end[0])/2
            mid_lon = (start[1] + end[1])/2
            midpoints.append((mid_lat, mid_lon))    

        # Plot the polygon with annotations
        # GeodesicService.plot_sector_with_annotations(polygon, point, midpoints, projections)

        if polygon.contains(point):
            return -min_distance    # Negative distance when inside the polygon
        return min_distance

    @staticmethod
    def create_polygon_from_coords(coords):
        """
        Create a Shapely polygon from a list of coordinates.
        Automatically converts coordinates from (lat, lon) to (lon, lat).
        """

        def sort_polygon_points(coords):
            """
            Sort the points of a polygon in counterclockwise order to ensure a regular polygon.
            """

            # Calculate the centroid
            coords_size = len(coords)
            centroid_lat = sum(p[0] for p in coords)/coords_size
            centroid_lon = sum(p[1] for p in coords)/coords_size

            # Sort points based on angle from the centroid
            def angle_from_centroid(point):
                lat, lon = point
                return math.atan2(lat - centroid_lat, lon - centroid_lon)

            sorted_coords = sorted(coords, key=angle_from_centroid)
            return sorted_coords
    
        sorted_coords = sort_polygon_points(coords)
        # Convert coordinates from (lat, lon) to (lon, lat)
        normalized_coords = [(coord[1], coord[0]) for coord in sorted_coords]
        normalized_coords.append(normalized_coords[0])  # Close the polygon
        sorted_coords.append(sorted_coords[0])  # Close the polygon

        return Polygon(normalized_coords), sorted_coords

    @staticmethod
    def f(line_coords, polygon_coords):
        """
        Use aircraft's current position to determine the horizontal slice of a sector. 
        """

        """curr_flight_level = current_position[-1] * 0.3048 # ft to meters

        print("Line coords initial: ", line_coords)
        line_coords_cartesian, polygon_coords_cartesian = [], [] 
        for p in line_coords:
            x, y, z = GeodesicService.xyz_from_lat_lon(p[0], p[1], curr_flight_level)
            line_coords_cartesian.append((x, y, z))
        print("Line coords cartesian: ", line_coords_cartesian)
        for p in polygon_coords:
            x, y, z = GeodesicService.xyz_from_lat_lon(p[0], p[1], curr_flight_level)
            polygon_coords_cartesian.append((x, y, z))
        print("Polygon coords cartesian: ", polygon_coords_cartesian)

        line_coords_normalized = [(lon, lat) for (lat, lon) in line_coords]
        linestring_latlon = LineString(line_coords_normalized)
        line = [(y, x, z) for (x, y, z) in line_coords_cartesian]
        linestring = LineString(line)
        print("Linestring: ", linestring)
        polygon_cartesian = [(x, y, z) for (x, y, z) in polygon_coords_cartesian]
        polygon_cartesian, _ = GeodesicService.create_polygon_from_coords(polygon_cartesian, curr_flight_level)
        print("Polygon cartesian: ", polygon_cartesian)
       
        # Calculate intersection
        intersection_cartesian = list(linestring.intersection(polygon_cartesian).coords)
        print("Intersection cartesian: ", intersection_cartesian)

        intersection_latlon = []
        for p in intersection_cartesian:
            lat, lon, height = GeodesicService.xyz_to_lat_lon(p[0], p[1], curr_flight_level)
            intersection_latlon.append((lat, lon, height))
        print("Intersection in latlon: ", intersection_latlon)

        polygon_latlon, _ = GeodesicService.create_polygon_from_coords(polygon_coords, curr_flight_level)
        # Plot the polygon
        poly_x, poly_y = zip(*polygon_latlon.exterior.coords)
        plt.plot(poly_x, poly_y, label="Polygon", color="blue")
        # Plot the line
        line_x, line_y = zip(*linestring_latlon.coords)
        plt.plot(line_x, line_y, label="Line", color="green")
        # Add labels and legend
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title("Line-Polygon Intersection")
        plt.legend()
        plt.grid(True)

        # Plot the intersection points
        if intersection_latlon:
            intersection_x, intersection_y = zip(*intersection_latlon)
            plt.scatter(intersection_y, intersection_x, color="red", label="Intersection Points")
            for i, point in enumerate(intersection_latlon):
                plt.text(point[0], point[1], str(i), color='black', fontsize=12, ha='center', va='center')
        
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()"""

        polygon, _ = GeodesicService.create_polygon_from_coords(polygon_coords)
        normalized_line_coords = [(coord[1], coord[0]) for coord in line_coords]
        line = LineString(normalized_line_coords) # TODO: create a method for line creation since it is used often everywhere
        print(f"polygon.area ={polygon.area}")
        """ # Visualization
        # Plot the polygon
        poly_x, poly_y = zip(*polygon.exterior.coords)
        plt.plot(poly_x, poly_y, label="Polygon", color="blue")
        # Plot the line
        line_x, line_y = zip(*line.coords)
        plt.plot(line_x, line_y, label="Line", color="green")
        # Add labels and legend
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title("Line-Polygon Intersection")
        plt.legend()
        plt.grid(True) """
        
        # Transforming the coordinates into a projected coordinate system
        """  projection = pyproj.Transformer.from_proj(
            pyproj.Proj("EPSG:4326"), 
            pyproj.Proj("EPSG:25832"),
            always_xy=True
        ) """
        
        projection = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32617", always_xy=True)
        inverse_projection = pyproj.Transformer.from_crs("EPSG:32617", "EPSG:4326", always_xy=True)

        projected_coords = [projection.transform(lon, lat) for lon, lat in polygon.exterior.coords]
        polygon = Polygon(projected_coords)

        #polygon = transform(projection.transform, polygon)
        line = transform(projection.transform, line)

        # Extracting the intersection
        intersection = line.intersection(polygon)

        # Inverse transformation
        """ inverse_projection = pyproj.Transformer.from_proj(
            pyproj.Proj("EPSG:25832"),
            pyproj.Proj("EPSG:4326"), 
            always_xy=True
        ) """
        intersection = transform(inverse_projection.transform, intersection)
        #print("Intersection points: ", intersection)
        list_of_intersection_candidates = []
        if intersection.geom_type == 'MultiLineString':
            for line in intersection.geoms:
                list_of_intersection_candidates.extend(list(line.coords))  # extend umjesto append
        elif intersection.geom_type == 'LineString':
            list_of_intersection_candidates = list(intersection.coords)  # bez dodatnih zagrada
            
        print(f"Final result structure: {list_of_intersection_candidates}")
        print(f"First element type: {type(list_of_intersection_candidates[0]) if list_of_intersection_candidates else 'Empty'}")
        return list_of_intersection_candidates

    @staticmethod
    def find_intersection_points(trajectory, polygon_coords):
        """
        Finds all intersection points of trajectory and a polygon.
        This method returns same values as previous project's Java methods for finding intersection points.
        """
        intersection_points_list = []
        
        for i in range(1, len(trajectory)):
            # Extract the next and previous trajectory points (lat, lon, time, height)
            next_lat, next_lon, _, _ = trajectory[i]
            prev_lat, prev_lon, _, _ = trajectory[i - 1]

            # Create a line segment between the next and previous trajectory points
            line_coords = [(prev_lat, prev_lon), (next_lat, next_lon)]

            # Check for intersection
            intersection_points = GeodesicService.line_sector_intersection(line_coords, polygon_coords)

            if intersection_points:
                # Add intersection points to the list, avoiding duplicates
                for point in intersection_points:
                    if point not in intersection_points_list:
                        intersection_points_list.append(point)
        
        return intersection_points_list

    @staticmethod
    def line_sector_intersection(line_coords, polygon_coords):
        """
        Find intersection points between a line segment and a polygon.
        Automatically converts coordinates to (lon, lat) if necessary.
        """

        # TODO: update with z coordinate (height)
        # TODO: could be replaced with .intersects()? This returns true/false.

        # Convert line coordinates from (lat, lon) to (lon, lat)
        normalized_line_coords = [(coord[1], coord[0]) for coord in line_coords]
        line = LineString(normalized_line_coords)
        polygon, _ = GeodesicService.create_polygon_from_coords(polygon_coords)
        intersection_candidates = line.intersection(polygon)
        print("Intersection candidates:", intersection_candidates)

        # Check if intersection_candidates is not empty
        if intersection_candidates.is_empty:
            list_of_intersection_candidates = []
        else:
            if intersection_candidates.geom_type == 'MultiLineString':
                list_of_intersection_candidates = [list(line.coords) for line in intersection_candidates.geoms]
        # If intersection_candidates is a LineString
            elif intersection_candidates.geom_type == 'LineString':
                list_of_intersection_candidates = [list(intersection_candidates.coords)]
        # Print the resulting list of intersection candidates
        print("List of intersection candidates:", list_of_intersection_candidates)
        intersection_candidates_expanded = []

        for inner_list in list_of_intersection_candidates:
            for candidate in inner_list:  # Iterating over each coordinate pair
                if isinstance(candidate, tuple):  # If candidate is a tuple, it's a coordinate pair
                    expanded_coords = list(candidate)  # Convert tuple to list
                    intersection_candidates_expanded.append(expanded_coords)
                else:
                    # If candidate is not a tuple, add it directly (this case shouldn't occur with valid input)
                    intersection_candidates_expanded.append([candidate])
        #print("Intersection candidates:", intersection_candidates)
        
        # For each intersection point, calculate the wgs84 distance to the projection of that point on the polygon edge.
        # Filter points that are around the edge by some tolerance.
        for point in intersection_candidates_expanded:
            print("\n")
            print("Point: ", Point(point))
            print("Distance to polygon: ", GeodesicService.calculate_distance_to_sector(point[1], point[0], polygon_coords))
            print("Polygon contains point: ", polygon.contains(Point(point)))
            print("Polygon contains point v2: ", GeodesicService.check_if_point_in_polygon(point[1], point[0], polygon_coords))
            print("Polygon contains point properly: ", polygon.contains_properly(Point(point)))
            print("Exterior contains point: ", polygon.exterior.contains(Point(point)))
            print("Exterior contains point properly: ", polygon.exterior.contains_properly(Point(point)))

        intersection = [
            (point[1], point[0]) for point in intersection_candidates_expanded if 
            -100 <= GeodesicService.calculate_distance_to_sector(point[1], point[0], polygon_coords) < 100 # 100 meters tolerance
        ]
        print("Intersection:", intersection)

        # Visualization
        # Plot the polygon
        poly_x, poly_y = zip(*polygon.exterior.coords)
        plt.plot(poly_x, poly_y, label="Polygon", color="blue")
        # Plot the line
        line_x, line_y = zip(*line.coords)
        plt.plot(line_x, line_y, label="Line", color="green")
        # Add labels and legend
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title("Line-Polygon Intersection")
        plt.legend()
        plt.grid(True)

        # Plot the intersection points
        if intersection:
            intersection_x, intersection_y = zip(*intersection)
            plt.scatter(intersection_y, intersection_x, color="red", label="Intersection Points")
            for i, point in enumerate(intersection):
                plt.text(point[0], point[1], str(i), color='black', fontsize=12, ha='center', va='center')
        
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
        
        return intersection
    
    @staticmethod
    def trajectory_sector_intersection(trajectory, polygon_coords):
        """
        Find all intersection points between a trajectory and a polygon.
        
        Args:
            trajectory (list of tuple): List of trajectory points in [(lat, lon, time, height), ...] format.
            polygon_coords (list of tuple): List of (lat, lon) tuples defining the polygon.

        Returns:
            list of tuples: List of (lat, lon) intersection points.
        """

        # TODO: update with z coordinate (height)

        # Create a polygon from the coordinates
        polygon, _ = GeodesicService.create_polygon_from_coords(polygon_coords)

        plt.clf()
        plt.close()

        
        # Plot the trajectory
        traj_lons = [point[1] for point in trajectory]
        traj_lats = [point[0] for point in trajectory]
        plt.plot(traj_lons, traj_lats, label="Trajectory", color="green", linestyle='--')

        intersections = GeodesicService.find_intersection_points(trajectory,polygon_coords)
        if intersections:
            intersection_y, intersection_x = zip(*intersections)

            # Plot intersection points directly on the polygon
            for point in intersections:
                y = point[1]
                x = point[0]
                plt.scatter(y, x, color="blue", zorder=5)  # Plot points on the polygon
                plt.text(y, x, f"({x:.2f}, {y:.2f})", fontsize=9, ha='right', va='bottom', color="blue")

        print("Intersections: ", intersections)
        # Add labels and legend
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title("Trajectory and Polygon Intersection")
        plt.grid(True)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.legend()
        plt.show()
        return intersections
    
    @staticmethod
    def find_exit_point(trajectory, polygon_coords, current_position):
        """
        Find the exit point(s) of a trajectory from a polygon.
        
        Args:
            trajectory (list of tuple): List of trajectory points in [(lat, lon, time, height), ...] format.
            polygon_coords (list of tuple): List of (lat, lon) tuples defining the polygon.
            current_position (tuple): Current position as (lat, lon).

        Returns:
            str: Exit point as a tuple (lat, lon).
        """

        # TODO: update with z coordinate (height)

        # Add current position to front, to create the first line segment
        trajectory.insert(0, (current_position[0], current_position[1], 0, 0))  # time and height are irrelevant
        #TODO: put this in separate method because the intersection points search is done in more methods
        with contextlib.redirect_stdout(io.StringIO()):
            exit_point_candidates = GeodesicService.find_intersection_points(trajectory,polygon_coords)
        
        # Plot the polygon - testing
        polygon, _ = GeodesicService.create_polygon_from_coords(polygon_coords)
        polygon_x, polygon_y = zip(*polygon.exterior.coords)
        plt.plot(polygon_x, polygon_y, label="Polygon", color="blue")

        # Plot current position
        plt.scatter(current_position[1], current_position[0], color="yellow", label="Current position", zorder=5)

        # Plot the trajectory (line segments)
        traj_lons = [point[1] for point in trajectory]
        traj_lats = [point[0] for point in trajectory]
        plt.plot(traj_lons, traj_lats, label="Trajectory", color="green", linestyle='--')

        # Add labels
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title("Trajectory and Polygon with Exit Point")
        plt.grid(True)
        plt.gca().set_aspect('equal', adjustable='box')

        if not exit_point_candidates:
            print("Exit point(s): ", exit_point_candidates)
            # Show the plot
            plt.legend()
            plt.show()
            return None
        
        if len(exit_point_candidates) == 1:
            last_lat, last_lon, _, _ = trajectory[-1]
            last_point = Point(last_lon, last_lat)

            if polygon.contains(last_point):
                print("Exit point(s): ", [])
                # Show the plot
                plt.legend()
                plt.show()
                return None
            else:
                print("Exit point(s): ", [(point[1], point[0]) for point in exit_point_candidates])
                # Plot the exit point
                plt.scatter(exit_point_candidates[0][1], exit_point_candidates[0][0], color="red", label="Exit Point", zorder=5)
                # Show the plot
                plt.legend()
                plt.show()
                return [(point[1], point[0]) for point in exit_point_candidates]    # (lat, lon) format
            
        if len(exit_point_candidates) > 1:  # TODO: If many exit points, check if a/c leaves and entries the sector in small time interval, this should be treated as not leaving the sector at all (i.e. we care only for the second exit point)
            first_lat, first_lon, _, _ = trajectory[0]
            first_point = Point(first_lon, first_lat)

            if polygon.covers(first_point):
                exit_points = exit_point_candidates[0::2] # Only odd intersection points are exit points
                print("Exit point(s): ", [(point[1], point[0]) for point in exit_points])
                for exit_point in exit_points:
                    # Plot the exit point
                    plt.scatter( exit_point[1],exit_point[0], color="red", label="Exit Point", zorder=5)
                # Show the plot
                plt.legend()
                plt.show()
                return [(point[1], point[0]) for point in exit_points]  # (lat, lon) format
            else:
                exit_points = exit_point_candidates[1::2] # Only even intersection points are exit points
                print("Exit point(s): ", [(point[1], point[0]) for point in exit_points])
                for exit_point in exit_points:
                    # Plot the exit point
                    plt.scatter( exit_point[1],exit_point[0], color="red", label="Exit Point", zorder=5)
                # Show the plot
                plt.legend()
                plt.show()
                return [(point[1], point[0]) for point in exit_points]  # (lat, lon) format
 
    @staticmethod
    def find_entry_point(trajectory, polygon_coords, current_position):
            
        """
        Find the entry point(s) of a trajectory into a polygon.

        Args:
            trajectory (list of tuple): List of trajectory points in [(lat, lon, time, height), ...] format.
            polygon_coords (list of tuple): List of (lat, lon) tuples defining the polygon.
            current_position (tuple): Current position as (lat, lon).

        Returns:
            tuple: Entry point as a (lat, lon) tuple.
        """

        # TODO: update with z coordinate (height)

        # Create a polygon from the coordinates
        polygon, _ = GeodesicService.create_polygon_from_coords(polygon_coords)
        
        trajectory.insert(0, (current_position[0], current_position[1], 0, 0))  # time and height are irrelevant
        # Iterate through the trajectory to find the first segment that intersects the polygon
        entry_point_candidates = GeodesicService.find_intersection_points(trajectory,polygon_coords)
        # Add current position to front, to create the first line segment
        current_point = Point(current_position[1], current_position[0])  # Convert to (lon, lat)

        # Clear any previous plots
        plt.clf()
        plt.close()

        # Plot the polygon
        polygon_x, polygon_y = zip(*polygon.exterior.coords)
        plt.plot(polygon_x, polygon_y, label="Polygon", color="blue")

        # Plot current position
        plt.scatter(current_position[1], current_position[0], color="yellow", label="Current position", zorder=5)

        # Plot the trajectory
        traj_lons = [point[1] for point in trajectory]
        traj_lats = [point[0] for point in trajectory]
        plt.plot(traj_lons, traj_lats, label="Trajectory", color="green", linestyle='--')

        # Add labels
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title("Trajectory and Polygon with Entry Point")
        plt.grid(True)
        plt.gca().set_aspect('equal', adjustable='box')

        # Preconditions: Current position must be outside the polygon
        if polygon.covers(current_point):
            print("Current position is within the sector, no entry point exists.")
            plt.legend()
            plt.show()
            return None
        # If no entry points found
        if not entry_point_candidates:
            print("No entry points found.")
            plt.legend()
            plt.show()
            return None

        
        lat, lon = entry_point_candidates[0][0], entry_point_candidates[0][1]  # Swap if necessary
        plt.scatter(lon, lat, color="red", label="Entry Point", zorder=5)

        # Show the plot
        plt.legend()
        plt.show()

        # Return the first entry point (in (lat, lon) format)
        return entry_point_candidates[0][0], entry_point_candidates[0][1]
    @staticmethod
    def calculate_time_to_point(lat1, lon1, lat2, lon2, speed):
        """
        Calculate the time in seconds required to travel between two points at a given speed.
        """

        distance = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)['s12']
        speed_m_per_s = speed * GeodesicService.c  # Convert knots to m/s
        print("Time to point [s]:", distance/speed_m_per_s)
        return distance/speed_m_per_s

    @staticmethod
    def distance_to_exit_point(trajectory, polygon_coords, current_position):
        """
        Calculate distance from current position to exit point
        """
        exit_point = GeodesicService.find_exit_point(trajectory, polygon_coords, current_position)[0]
        exit_point_data = {}
        print(f"exit point je: {exit_point}")  
        exit_point_data["lat"],exit_point_data["lon"] =  current_position[1], current_position[0]
        distance = Geodesic.WGS84.Inverse(current_position[1], current_position[0], exit_point[1], exit_point[0])['s12']
        exit_point_data["distance"] = distance
        return exit_point_data
    
    @staticmethod        
    def military_sector_intersection(trajectory, current_position, coords,direct_to_point):
        """
        Calculates the value of flag variable used in Task 8.1 for detecting direct-to points
        """
        if direct_to_point == (0, 0):
            direct_to_lat,direct_to_lon = GeodesicService.find_exit_point(trajectory, coords, current_position)[0]
        else:    
            direct_to_lat,direct_to_lon = direct_to_point[1],direct_to_point[0] 
            
        current_lon,current_lat  = current_position[1], current_position[0]
        azi = Geodesic.WGS84.Inverse(current_lat, current_lon, direct_to_lat, direct_to_lon)["azi1"]
        flag = 1
        for i in range(0, len(coords), 1):
            j = i + 1
            if i == len(coords) - 1:
                j = 0
            lat1, lon1 = coords[i]
            lat2, lon2 = coords[j]
            azi1 =  Geodesic.WGS84.Inverse(current_lat, current_lon, lat1, lon1)["azi1"]
            azi2 =  Geodesic.WGS84.Inverse(current_lat, current_lon, lat2, lon2)["azi1"]
			
            if (azi1 < 0 and azi2 > 0) or (azi1 > 0 and azi2 < 0):
                if abs(azi1) < 90 and abs(azi2) < 90:
                    if azi1 < azi2:
                        if azi > azi1 and azi < azi2:
                            flag = 0
                    else:
                        if azi < azi1 and azi > azi2:
                            flag = 0
                elif (abs(azi1) > 90 and abs(azi2) > 90) or (abs(azi1) < 90 and abs(azi2) > 90) or (abs(azi1) > 90 and abs(azi2) < 90):
                    if azi1 < azi2:
                        if azi < azi1 and azi > azi2:
                            flag = 0
                    else:
                        if azi > azi1 and azi < azi2:
                            flag = 0
            else:
                if azi1 < azi2:
                    if azi > azi1 and azi < azi2:
                        flag = 0
                else:
                    if azi < azi1 and azi > azi2:
                        flag = 0
                        
        return flag
    
    @staticmethod
    def flying_towards_exit_point(heading, trajectory, polygon_coords, current_position):
        """
        Checks if the aircraft is flying towards the exit point.
        
        Parameters:
        ----------
        heading : float
            The current heading of the aircraft, received from flight data.
        trajectory : list
            The flight trajectory, typically a list of waypoints.
        polygon_coords : list
            Coordinates defining the military sector or restricted airspace.
        current_position : tuple
            The current position of the aircraft as (latitude, longitude).

        Returns:
        -------
        1/0
            True if the aircraft is flying towards the exit point, False otherwise.
        """
        
        current_lon,current_lat  = current_position[1], current_position[0]
        exit_point_data = GeodesicService.distance_to_exit_point(trajectory, polygon_coords, current_position)
        exit_point_lat,exit_point_lon = exit_point_data["lat"],exit_point_data["lon"] 
        azi1 = Geodesic.WGS84.Inverse(current_lat, current_lon, exit_point_lat, exit_point_lon)["azi1"]
        radius = 2.5 * 1825  # meters
        distance_to_exit_point = exit_point_data["distance"]
        alpha_in_pi = math.atan(radius / distance_to_exit_point)
        alpha_in_deg = math.degrees(alpha_in_pi)  # Convert radians to degrees
        if azi1 + alpha_in_deg >= heading >= azi1 - alpha_in_deg:
            return 1 #"Aircraft is flying towards the exit point."
        else:
            return 0 #"Aircraft is not flying towards the exit point."
        
    @staticmethod
    def route_point_closest_to_exit_point(trajectory, polygon_coords, current_position):
        """
        Finds the route point closest to the exit point and returns distances for all points.
        
        Args:
            trajectory: List of tuples (lat, lon) - route coordinates
            polygon_coords: Polygon coordinates
            current_position: Current position (lat, lon)
            
        Returns:
            tuple: (closest_point, OrderedDict with distances for all points)
        """
        # Find the exit point
        trajectory_copy = list(trajectory)
        with contextlib.redirect_stdout(io.StringIO()):
            exit_points = GeodesicService.find_exit_point(trajectory_copy, polygon_coords, current_position)
        if not exit_points:
            raise ValueError("No exit points found")
        #print(exit_points)
        exit_lat, exit_lon = exit_points[-1][1], exit_points[-1][0]
        distances = OrderedDict()

        # Iterate through all route segments
        for i in range(1, len(trajectory)):
            prev_point = trajectory[i-1]
            next_point = trajectory[i]
            
            # Extract coordinates (assumption: lat is the first element, lon is the second)
            prev_lat, prev_lon = prev_point[0], prev_point[1]
            next_lat, next_lon = next_point[0], next_point[1]
            
            # Create a line segment
            segment = [(prev_lat, prev_lon), (next_lat, next_lon)]
            
            # Check intersection with the sector
            with contextlib.redirect_stdout(io.StringIO()):
                if GeodesicService.line_sector_intersection(segment, polygon_coords):
                    # Add both segment points for consideration
                    for point in (prev_point, next_point):
                        point_lat, point_lon = point[0], point[1]
                        print(exit_lat, exit_lon, point_lat, point_lon)
                        distance = GeodesicService.geodesic_distance(
                            exit_lat, exit_lon,
                            point_lat, point_lon
                        )
                        distances[point] = distance
        #print(exit_lat, exit_lon, point_lat, point_lon)
        if not distances:
            raise ValueError("No segments intersect the sector")
        
        closest_point = min(distances, key=distances.get)
        return closest_point, distances

    @staticmethod
    def route_point_closest_to_entry_point(trajectory, polygon_coords, current_position):
        """
        Finds the route point closest to the entry point and returns distances for all points.
        
        Args:
            trajectory: List of tuples (lat, lon) - route coordinates
            polygon_coords: Polygon coordinates
            current_position: Current position (lat, lon)
            
        Returns:
            tuple: (closest_point, OrderedDict with distances for all points)
        """
        # Find the entry point
        #trajectory_copy = list(trajectory)
        with contextlib.redirect_stdout(io.StringIO()):
            entry_points = GeodesicService.find_entry_point(trajectory, polygon_coords, current_position)
        if not entry_points:
            raise ValueError("No entry points found")
        print(f"entry: {entry_points}")
        entry_lat, entry_lon = entry_points[0], entry_points[1]
        distances = OrderedDict()

        # Iterate through all route segments
        for i in range(1, len(trajectory)):
            
            prev_point = trajectory[i-1]
            next_point = trajectory[i]
            
            # Extract coordinates (assumption: lat is the first element, lon is the second)
            prev_lat, prev_lon = prev_point[0], prev_point[1]
            next_lat, next_lon = next_point[0], next_point[1]
            
            # Create a line segment
            segment = [(prev_lat, prev_lon), (next_lat, next_lon)]
            
            # Check intersection with the sector
            with contextlib.redirect_stdout(io.StringIO()):
                intersection = GeodesicService.line_sector_intersection(segment, polygon_coords)

            # The rest of the code prints normally since we are not inside the redirect block
            if intersection:
                for point in (prev_point, next_point):
                    
                    if point == (*current_position,0,0):
                        continue
                    point_lat, point_lon = point[0], point[1]
                    #print(entry_lat, entry_lon, point_lat, point_lon)
                    distance = GeodesicService.geodesic_distance(
                        entry_lat, entry_lon,
                        point_lat, point_lon
                    )
                    distances[point] = distance
                    #print(entry_lat, entry_lon, point_lat, point_lon)
                    if not distances:
                        raise ValueError("No segments intersect the sector")
        
        closest_point = min(distances, key=distances.get)
        return closest_point, distances