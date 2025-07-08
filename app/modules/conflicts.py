class ConflictDetection:
    def __init__(self):
        self.detections = {} # key: timestamp, value: list of detected conflicts

    def _extract_repo_cd_data(self, data_repository, timestamp):
        """
        Extracts data from the data repository for conflict detection (CD).
        Gathers current position, flight level, predicted trajectory, cruise level,
        waypoints, speed, and airport information for each flight at the given timestamp.
        """
        state_data = data_repository.get(timestamp, {})
        cd_data = {timestamp: []}

        for flight_key, flight_data in state_data.items():
            rtg = flight_data["Flight"]["routeTrajectoryGroup"]

            # Skip if 'predicted' branch doesn't exist
            if "predicted" not in rtg:
                continue

            current_point = rtg.get("current", {}).get("element", {}).get("point4D", {})
            current_pos = current_point.get("position", {}).get("pos")
            current_fl = current_point.get("level", {}).get("flightLevel")

            # Predicted trajectory points
            predicted_points = rtg["predicted"].get("points", [])

            cleared_branch = rtg.get("agreed", [])
            cleared_waypoints = [entry.get("element", {}).get("point4D", {}) for entry in cleared_branch]
            # Find the dict in the list that contains the "routeInformation" key
            cleared_route_info_dict = next(
                (item for item in cleared_branch if "routeInformation" in item),
                {}
            )
            # Access flight level from the route information in cleared branch
            cruise_level_cleared = cleared_route_info_dict.get("routeInformation", {}).get("cruisingLevel", {}).get("flightLevel")

            desired_branch = rtg.get("desired", [])
            # Find the dict in the list that contains the "routeInformation" key
            desired_route_info_dict = next(
                (item for item in desired_branch if "routeInformation" in item),
                {}
            )
            # Access flight level from the route information in desired branch
            cruise_level_desired = desired_route_info_dict.get("routeInformation", {}).get("cruisingLevel", {}).get("flightLevel")

            flight = flight_data["Flight"]
            exit_point = flight.get("enRoute", {}).get("exitPoint")

            departure_location_indicator = flight.get("departure", {}).get("departureAerodrome", {}).get("locationIndicator")
            arrival = flight.get("arrival", {})
            arrival_location_indicator = arrival.get("destinationAerodrome", {}).get("locationIndicator")

            alt_arrival = arrival.get("alternateAerodromeAlternate", {})
            alt_arrival_location_indicator = alt_arrival.get("locationIndicator") if alt_arrival else None

            cd_data[timestamp].append({
                "callsign": flight_key,
                "current_pos": current_pos,
                "current_fl": current_fl,
                "predicted_points": predicted_points,
                "cruise_level_cleared": cruise_level_cleared,
                "cruise_level_desired": cruise_level_desired,
                "exit_point": exit_point,
                "cleared_waypoints": cleared_waypoints,
                "departure_aerodrome_indicator": departure_location_indicator,
                "arrival_aerodrome_indicator": arrival_location_indicator,
                "arrival_aerodrome_alternate_indicator": alt_arrival_location_indicator
            })

        return cd_data

    
    def _find_conflicts(self, records):
        """
        This method implements the logic to identify conflicts based on the data from repository.
        """
        pass
    
    def detect(self, data_repository, last_timestamp):
        """
        Detect conflicts in the data repository based on the last timestamp.
        This method should implement the logic to identify conflicts in the data.
        """
        # Example conflict detection logic (to be replaced with actual implementation)
        for timestamp, records in data_repository.items():
            if timestamp > last_timestamp:
                # Check for conflicts in records
                conflicts = self._find_conflicts(records)
                if conflicts:
                    self.detections[timestamp] = conflicts


class ConflictResolution:
    def __init__(self):
        self.resolutions = {} # key: timestamp, value: list of resolutions for conflicts

    def resolve(self, conflict_detections):
        """
        Resolve conflicts based on the detected conflicts.
        This method implement the logic to resolve conflicts.
        """
        pass