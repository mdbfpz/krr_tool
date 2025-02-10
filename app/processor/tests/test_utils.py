import unittest
import sys
import os

# Add the project root to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from app.processor.utils import GeodesicService


class TestGeodesic(unittest.TestCase):

    """def test_check_if_point_in_circle(self):
        # Eiffel Tower within 5 km of Louvre

        self.assertTrue(GeodesicService.check_if_point_in_circle(48.858844, 2.294351, 48.860611, 2.337644, 5000))
    
    def test_calculate_tolerance_azi(self):
        # Tolerance for a 1 km radius

        tolerance = GeodesicService.calculate_tolerance_azi(48.858844, 2.294351, 48.860611, 2.337644, 1000)
        self.assertGreater(tolerance, 0)
    
    def test_xyz_from_lat_lon(self):
        # Test the conversion from geographic to Cartesian coordinates.

        # Mount Everest's summit (lat, lon in decimal degrees, height in meters)
        lat, lon, height = 27.9881, 86.9250, 8848

        # Convert to Cartesian coordinates
        x, y, z = GeodesicService.xyz_from_lat_lon(lat, lon, height)

        # Expected Cartesian coordinates (values are approximate)
        self.assertAlmostEqual(x, 302769.89, delta=2)
        self.assertAlmostEqual(y, 5636025.47, delta=2)
        self.assertAlmostEqual(z, 2979493.09, delta=1)

    def test_xyz_to_lat_lon(self):
        # Test the conversion from Cartesian to geographic coordinates.

        # Cartesian coordinates of Mount Everest's summit
        x, y, z = 302769.8, 5636027.0, 2979493.087

        # Convert to geographic coordinates (lat, lon in decimal degrees, height in meters)
        lat, lon, height = GeodesicService.xyz_to_lat_lon(x, y, z)

        # Expected geographic coordinates (latitude and longitude in decimal degrees)
        self.assertAlmostEqual(lat, 27.9881, delta=0.0001)
        self.assertAlmostEqual(lon, 86.9250, delta=0.0001)
        self.assertAlmostEqual(height, 8848, delta=2)
    
    def test_check_if_point_in_polygon(self):
        # Test Case 1: Point inside the polygon
        polygon_coords = [(40.6892, -74.0455), (40.7003, -74.0204), (40.6793, -74.0455)]
        self.assertTrue(GeodesicService.check_if_point_in_polygon(40.687, -74.04, polygon_coords))

        # Test Case 2: Point outside the polygon
        polygon_coords = [(40.6892, -74.0455), (40.7003, -74.0204), (40.6793, -74.0455)]
        self.assertFalse(GeodesicService.check_if_point_in_polygon(40.710, -74.05, polygon_coords))

        # Test Case 3: Point exactly on the boundary of the polygon
        polygon_coords = [(40.6892, -74.0455), (40.7003, -74.0204), (40.6793, -74.0455)]
        self.assertFalse(GeodesicService.check_if_point_in_polygon(40.6892, -74.0455, polygon_coords))  # On the boundary

        # Test Case 4: Point inside a different polygon (larger)
        polygon_coords = [(40.7300, -73.9350), (40.7400, -73.9250), (40.7250, -73.9150), (40.7150, -73.9250)]
        self.assertTrue(GeodesicService.check_if_point_in_polygon(40.7300, -73.9250, polygon_coords))  # Inside the polygon

        # Test Case 5: Point far outside the polygon
        polygon_coords = [(40.7300, -73.9350), (40.7400, -73.9250), (40.7250, -73.9150), (40.7150, -73.9250)]
        self.assertFalse(GeodesicService.check_if_point_in_polygon(41.0, -74.0, polygon_coords))  # Far outside
    
    def test_calculate_distance_to_sector(self):
        # Test Case 1: Point inside the polygon
        polygon_coords = [(39.967778, -83.024722), (39.949722, -75.164167), (38.916667, -77.070556), (41.464167, -81.664444)]
        distance = GeodesicService.calculate_distance_to_sector(40.2, -79.1, polygon_coords)
        self.assertLess(distance, 0)

        # Test Case 2: Point outside the polygon (close to the boundary)
        polygon_coords = [(39.967778, -83.024722), (39.949722, -75.164167), (38.916667, -77.070556), (41.464167, -81.664444)]
        distance = GeodesicService.calculate_distance_to_sector(40.720, -74.05, polygon_coords)
        self.assertGreater(distance, 0)

        # Test Case 3: Point exactly on the boundary of the polygon
        polygon_coords = [(39.967778, -83.024722), (39.949722, -75.164167), (38.916667, -77.070556), (41.464167, -81.664444)]
        distance = GeodesicService.calculate_distance_to_sector(39.967778, -83.024722, polygon_coords)
        self.assertEqual(distance, 0)

        # Test Case 4: Point far outside the polygon
        polygon_coords = [(39.967778, -83.024722), (39.949722, -75.164167), (38.916667, -77.070556), (41.464167, -81.664444)]
        distance = GeodesicService.calculate_distance_to_sector(35.0, -85.0, polygon_coords)
        self.assertGreater(distance, 0)

        # Test Case 5: Point very close to one of the polygon's vertices
        polygon_coords = [(39.967778, -83.024722), (39.949722, -75.164167), (38.916667, -77.070556), (41.464167, -81.664444)]
        distance = GeodesicService.calculate_distance_to_sector(39.91, -83.02, polygon_coords)
        self.assertGreater(distance, 0)
    
    def test_create_polygon_from_coords(self):
        # Create a polygon

        polygon, _ = GeodesicService.create_polygon_from_coords([(39.967778, -83.024722), (39.949722, -75.164167), (38.916667, -77.070556), (41.464167, -81.664444)])
        self.assertTrue(polygon.is_valid)
    """
    def test_f(self):
        polygon_coords = [
            (39.9726058, -83.0202268),
            (39.9583192, -75.1627617),
            (38.9216105, -77.0526693),
            (41.4483224, -81.6853763)
        ]

        line_coords = [(39.3983302, -82.6424770), (40.7001176, -74.0349471)]
        current_position = (0, 0, 1000)
        GeodesicService.f(line_coords, polygon_coords, current_position)
        d = GeodesicService.geodesic_distance(lat1=40.445643882658466, lon1=-76.87341736600898, lat2=40.3474556, lon2=-76.8670986)
        print(d)
    """
        GeodesicService.f(line_coords, polygon_coords)
        
    
    def test_line_sector_intersection(self):
        polygon_coords = [#konkavan
            (37.9726058, -83.0202268),
            (40.5583192, -75.00),
            (39.9216105, -79.0526693),
            (41.4483224, -81.6853763)
        ]
        polygon_coords2 = [
            (39.9726058, -83.0202268),
            (39.9583192, -75.1627617),
            (38.9216105, -77.0526693),
            (41.4483224, -81.6853763)
        ]
        

        # Test Case 1: Intersection between line and polygon
        line_coords = [(38.7983302, -82.6424770), (40.7001176, -75.50349471)]
        intersection = GeodesicService.line_sector_intersection(line_coords, polygon_coords)
        #self.assertEqual(len(intersection), 3) 

        # Test Case 2: Line completely outside polygon (no intersection)
        line_coords = [(42.0, -85.0), (43.0, -84.0)]
        intersection = GeodesicService.line_sector_intersection(line_coords, polygon_coords)
        #self.assertEqual(len(intersection), 0)

        # Test Case 3: Line coincides with two of the polygon edges
        line_coords = [(39.967778, -83.024722), (39.949722, -75.164167)]
        intersection = GeodesicService.line_sector_intersection(line_coords, polygon_coords)
        #self.assertEqual(len(intersection), 2)

        # Test Case 4: Line starts inside the polygon and exits
        line_coords = [(40.0, -78.0), (41.0, -73.0)]
        intersection = GeodesicService.line_sector_intersection(line_coords, polygon_coords)
       # self.assertEqual(len(intersection), 1)

        # Test Case 5: Line starts outside and ends on the polygon vertex
        line_coords = [(38.0, -80.0), (39.949722, -75.164167)]
        intersection = GeodesicService.line_sector_intersection(line_coords, polygon_coords)
        #self.assertEqual(len(intersection), 2)

        # Test Case 6: Line starts outside and ends inside the polygon
        line_coords = [(38.0, -80.0), (39.9, -77.0)]
        intersection = GeodesicService.line_sector_intersection(line_coords, polygon_coords)
        self.assertEqual(len(intersection), 1)
    
        #self.assertEqual(len(intersection), 1)
    """
    def test_trajectory_sector_intersection(self):
        polygon_coords = [#konkavan
            (37.9726058, -83.0202268),
            (40.5583192, -75.00),
            (39.9216105, -79.0526693),
            (41.4483224, -81.6853763)
        ]
        
        """ polygon_coords = [
        (39.967778, -83.024722), 
        (39.949722, -75.164167), 
        (38.916667, -77.070556), 
        (41.464167, -81.664444)
        ]     """
        # Case 1: Trajectory starts and ends entirely outside the sector, zero intersections
        trajectory = [
            (39.0, -83.5, 1200, 13000),
            (40.0, -84.0, 1205, 13000), 
            (42.0, -82.0, 1210, 13000),
            (40.5, -75.0, 1210, 13000),
            (39.0, -75.0, 1210, 13000)
        ]
        trajectoryK = [
            (38.58, -83.27, 1200, 13000),
            (40.53, -81.0, 1205, 13000), 
            (42.13, -81.50, 1210, 13000),
            (41.25, -76.83, 1210, 13000),
            (38.8, -76.30, 1210, 13000)
        ]
        #intersections = GeodesicService.trajectory_sector_intersection(trajectoryK, polygon_coords)
        """
        #self.assertEqual(len(intersections), 0)

        # Case 2: Trajectory starts outside the sector, ends inside the sector, one intersection
        trajectory = [
            (39.0, -83.5, 1200, 13000),
            (40.5, -81.0, 1205, 13000)
        ]
        intersections = GeodesicService.trajectory_sector_intersection(trajectory, polygon_coords)
        #self.assertEqual(len(intersections), 1)
        
        # Case 3: Trajectory starts inside the sector and exits, ends outside, one intersection
        trajectory = [
            (40.0, -81.0, 1200, 13000),
            (40.5, -75.0, 1205, 13000)
        ]
        intersections = GeodesicService.trajectory_sector_intersection(trajectory, polygon_coords)
        #self.assertEqual(len(intersections), 1)

        # Case 4: Trajectory starts outside, passes through, and ends outside the sector, two intersection
        trajectory = [
            (39.4, -82.6, 1200, 13000),
            (39.8, -80.035, 1205, 13000), 
            (40.0, -77.035, 1210, 13000), 
            (40.7, -74.035, 1215, 13000)
        ]
        intersections = GeodesicService.trajectory_sector_intersection(trajectory, polygon_coords)
        #self.assertEqual(len(intersections), 2)

        # Case 5: Trajectory starts outside, exits and enters the sector multiple times  
        trajectory = [
            (39.4, -82.0, 1200, 13000),
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1215, 13000)
        ]
        intersections = GeodesicService.trajectory_sector_intersection(trajectory, polygon_coords)
        #self.assertEqual(len(intersections), 3)

        # Case 6: Trajectory starts outside, exits and enters the sector multiple times  
        trajectory = [
            (39.5, -82.0, 1200, 13000),
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1215, 13000),
            (40.5, -76.0, 1215, 13000)
        ]
        intersections = GeodesicService.trajectory_sector_intersection(trajectory, polygon_coords)
        #self.assertEqual(len(intersections), 4)

        # Case 7: Trajectory starts inside, exits and enters the sector multiple times
        trajectory = [
            (40.0, -82.0, 1200, 13000),
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1215, 13000)
        ]
        intersections = GeodesicService.trajectory_sector_intersection(trajectory, polygon_coords)
        #self.assertEqual(len(intersections), 2)

        # Case 8: Trajectory starts inside, exits and enters the sector multiple times
        trajectory = [
            (40.0, -82.0, 1200, 13000),
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1215, 13000),
            (40.5, -76.0, 1215, 13000)
        ]
        intersections = GeodesicService.trajectory_sector_intersection(trajectory, polygon_coords)
        #self.assertEqual(len(intersections), 3)
    """
    def test_find_exit_point(self):
        polygon_coords = [#konkavan
            (37.9726058, -83.0202268),
            (40.5583192, -75.00),
            (39.9216105, -79.0526693),
            (41.4483224, -81.6853763)
        ]

        """ # Case 1: Trajectory starts and ends entirely outside the sector, zero intersections
        polygon_coords = [
            (39.967778, -83.024722), 
            (39.949722, -75.164167), 
            (38.916667, -77.070556), 
            (41.464167, -81.664444)
        ]
        trajectory = [
            (40.0, -84.0, 1205, 13000), 
            (42.0, -82.0, 1210, 13000),
            (40.5, -75.0, 1210, 13000),
            (39.0, -75.0, 1210, 13000)
        ]
        current_position = (39.0, -83.5)
        exit_point = GeodesicService.find_exit_point(trajectory, polygon_coords, current_position)
        self.assertIsNone(exit_point, "Exit point should be None for a trajectory entirely outside the polygon.")
        # Case 2: Trajectory starts outside the sector, ends inside the sector, one intersection
        trajectory = [
            (40.5, -81.0, 1205, 13000)
        ]
        current_position = (39.0, -83.5)
        exit_point = GeodesicService.find_exit_point(trajectory, polygon_coords, current_position)
        self.assertIsNone(exit_point, "Exit point should be None for a trajectory that only touches the polygon boundary.")
        
         # Case 3: Trajectory starts inside the sector and exits, ends outside, one intersection
        trajectory = [
            (40.5, -75.0, 1205, 13000)
        ]
        current_position = (40, -81.0)
        exit_points = GeodesicService.find_exit_point(trajectory, polygon_coords, current_position)
        self.assertEqual(len(exit_points), 1)
        
        # Case 4: Trajectory starts outside, passes through, and ends outside the sector, two intersection
        trajectory = [
            (39.8, -80.035, 1205, 13000), 
            (40.0, -77.035, 1210, 13000), 
            (40.7, -74.035, 1215, 13000)
        ]
        current_position = (39.4, -82.6)
        exit_points = GeodesicService.find_exit_point(trajectory, polygon_coords, current_position)
        self.assertEqual(len(exit_points), 1)
        # Case 5: Trajectory starts outside, exits and enters the sector multiple times  
        trajectory = [
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1215, 13000)
        ]
        current_position = (39.4, -82.0)
        exit_points = GeodesicService.find_exit_point(trajectory, polygon_coords, current_position)
        self.assertEqual(len(exit_points), 1)

        # Case 6: Trajectory starts outside, exits and enters the sector multiple times  
        trajectory = [
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1215, 13000),
            (40.5, -76.0, 1215, 13000)
        ]
        current_position = (39.5, -82.0)
        exit_points = GeodesicService.find_exit_point(trajectory, polygon_coords, current_position)
        self.assertEqual(len(exit_points), 2)

        # Case 7: Trajectory starts inside, exits and enters the sector multiple times
        trajectory = [
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1215, 13000)
        ]
        current_position = (40.0, -82.0)
        exit_points = GeodesicService.find_exit_point(trajectory, polygon_coords, current_position)
        self.assertEqual(len(exit_points), 1)

        # Case 8: Trajectory starts inside, exits and enters the sector multiple times
        trajectory = [
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (39.0, -79.0, 1215, 13000),
            (40.5, -76.0, 1215, 13000)
        ]
        current_position = (40.0, -82.0)
        exit_points = GeodesicService.find_exit_point(trajectory, polygon_coords, current_position)
        #self.assertEqual(len(exit_points), 2)
        """
        
        # Case 9: Trajectory goes underneath
        trajectory = [
            
            (38.5, -80.5, 1202, 12500),  # Move towards the polygon
            (38.8, -79.0, 1205, 13000),  # Enter the polygon
            (39.3, -78.3, 1210, 13000),  # Move within the polygon
            (39.8, -76.40, 1215, 13000),  # Exit the polygon
            (41.0, -75.75, 1220, 13500)   # Additional point
            
        ]
        current_position = (39.2, -83.3)
        #exit_points = GeodesicService.find_exit_point(trajectory, polygon_coords, current_position)

        
        # Case 9: Trajectory starts outside, overlaps with the segment, ends outside
        # TODO: implement this test case - Kika: a/c can't fly along the boundary, horizontal distance to the boundary must be >= 2.5 NM
    
    def test_find_entry_point(self):
        polygon_coords = [#konkavan
            (37.9726058, -83.0202268),
            (40.5583192, -75.00),
            (40.6, -76.8),
            (40.40, -76.50),
        ]
        """
        # Case 1: Trajectory starts and ends entirely outside the sector, zero intersections
        trajectory = [
            (40.0, -84.0, 1205, 13000), 
            (42.0, -82.0, 1210, 13000),
            (40.5, -75.0, 1210, 13000),
            (39.0, -75.0, 1215, 13000)
        ]
        current_position = (39.0, -83.5)
        entry_point = GeodesicService.find_entry_point(trajectory, polygon_coords, current_position)

        # Case 2: Trajectory starts outside the sector, ends inside the sector, one intersection
        trajectory = [
            (39.95, -83.0, 1200, 13000), 
            (39.95, -80.0, 1205, 13000)
        ]
        current_position = (39.0, -84.0)
        entry_point = GeodesicService.find_entry_point(trajectory, polygon_coords, current_position)
        
        # Case 3: Trajectory starts inside the sector and exits, ends outside, one intersection -> no entry point
        trajectory = [
            (39.95, -81.0, 1200, 13000), 
            (40.5, -75.0, 1205, 13000)
        ]
        current_position = (39.9, -81.0)
        entry_point = GeodesicService.find_entry_point(trajectory, polygon_coords, current_position)

        # Case 4: Trajectory starts outside, passes through, and ends outside the sector, two intersection
        trajectory = [
            (39.39, -82.64, 1200, 13000), 
            (39.8, -80.035, 1205, 13000), 
            (40.0, -77.035, 1210, 13000), 
            (40.7, -74.035, 1215, 13000)
        ]
        current_position = (39.0, -83.0)
        entry_point = GeodesicService.find_entry_point(trajectory, polygon_coords, current_position)

        # Case 5: Trajectory starts outside, exits and enters the sector multiple times  
        trajectory = [
            (39.4, -82.0, 1200, 13000), 
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1220, 13000)
        ]
        current_position = (39.0, -83.0)
        entry_point = GeodesicService.find_entry_point(trajectory, polygon_coords, current_position)
        
        # Case 6: Trajectory starts outside, exits and enters the sector multiple times  
        trajectory = [
            (39.4, -82.0, 1200, 13000), 
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1220, 13000),
            (40.5, -76.0, 1225, 13000)
        ]
        current_position = (39.0, -83.0)
        entry_point = GeodesicService.find_entry_point(trajectory, polygon_coords, current_position)
        
        # Case 7: Trajectory starts inside, exits and enters the sector multiple times -> no entry point
        trajectory = [
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1225, 13000)
        ]
        current_position = (40.0, -82.0)
        entry_point = GeodesicService.find_entry_point(trajectory, polygon_coords, current_position)
        

        # Case 8: Trajectory starts inside, exits and enters the sector multiple times -> no entry point
        trajectory = [
            (40.3, -81.0, 1205, 13000), 
            (41.5, -79.5, 1210, 13000), 
            (41.0, -77.0, 1215, 13000),
            (40.0, -79.0, 1220, 13000),
            (40.5, -76.0, 1225, 13000)
        ]
        current_position = (40.0, -82.0)
        entry_point = GeodesicService.find_entry_point(trajectory, polygon_coords, current_position)
        """
        # Case 9: Trajectory goes underneath
        trajectory = [
            
            (38.5, -80.5, 1202, 12500),  # Move towards the polygon
            (38.8, -79.0, 1205, 13000),  # Enter the polygon
            (39.3, -78.3, 1210, 13000),  # Move within the polygon
            (39.8, -76.40, 1215, 13000),  # Exit the polygon
            (41.0, -75.75, 1220, 13500)   # Additional point
            
        ]
        current_position = (39.2, -83.3)
        #entry_point = GeodesicService.find_entry_point(trajectory, polygon_coords, current_position)

    def test_calculate_time_to_point(self):

        time = GeodesicService.calculate_time_to_point(40.6413, -73.7781, 40.7769, -73.8740, 200)
        #self.assertGreater(time, 0)
    
    def test_distance_to_exit_point(self):
        polygon_coords = [#konkavan
            (37.9726058, -83.0202268),
            (40.5583192, -75.00),
            (40.6, -76.8),
            (40.40, -76.50),
        ]
        trajectory = [
            
            (38.5, -80.5, 1202, 12500),  # Move towards the polygon
            (38.8, -79.0, 1205, 13000),  # Enter the polygon
            (39.3, -78.3, 1210, 13000),  # Move within the polygon
            (39.8, -76.40, 1215, 13000),  # Exit the polygon
            (41.0, -75.75, 1220, 13500)   # Additional point
            
        ]
        current_position = (39.2, -83.3)
        #distance_to_exit_point = GeodesicService.distance_to_exit_point(trajectory, polygon_coords, current_position)
        #print(distance_to_exit_point)

    def test_military_sector_intersection(self):
        polygon_coords = [#konkavan
            (37.9726058, -83.0202268),
            (40.5583192, -75.00),
            (40.6, -76.8),
            (40.40, -76.50),
        ]
        trajectory = [
            
            (38.5, -80.5, 1202, 12500),  # Move towards the polygon
            (38.8, -79.0, 1205, 13000),  # Enter the polygon
            (39.3, -78.3, 1210, 13000),  # Move within the polygon
            (39.8, -76.40, 1215, 13000),  # Exit the polygon
            (41.0, -75.75, 1220, 13500)   # Additional point
            
        ]
        current_position = (39.2, -83.3)
        direct_to_point = (0, 0)
        #military_sector_intersection = GeodesicService.military_sector_intersection(trajectory, current_position, polygon_coords,direct_to_point)
        #print(military_sector_intersection)
        
    def test_flying_towards_exit_point(self):
        polygon_coords = [#konkavan
            (37.9726058, -83.0202268),
            (40.5583192, -75.00),
            (40.6, -76.8),
            (40.40, -76.50),
        ]
        trajectory = [
            
            (38.5, -80.5, 1202, 12500),  # Move towards the polygon
            (38.8, -79.0, 1205, 13000),  # Enter the polygon
            (39.3, -78.3, 1210, 13000),  # Move within the polygon
            (39.8, -76.40, 1215, 13000),  # Exit the polygon
            (41.0, -75.75, 1220, 13500)   # Additional point
            
        ]
        current_position = (39.2, -83.3)
        heading = 0.7        
        flying_towards_exit_point_indicator = GeodesicService.flying_towards_exit_point(heading,trajectory, polygon_coords, current_position)
        print(flying_towards_exit_point_indicator)
if __name__ == "__main__":
    unittest.main()
