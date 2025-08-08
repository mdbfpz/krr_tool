from datetime import datetime, timedelta

class ConflictDetection:
    def __init__(self):
        self.detections = {} # key: timestamp, value: list of detected conflicts

    def _find_most_recent_predicted_trajectory(self, data_repository, flight_key):
        """
        Finds the most recent predicted trajectory for a flight from the data repository.
        This method retrieves the latest predicted trajectory points for a given flight key.
        """
        # Get all timestamps from the data repository
        timestamps = sorted(data_repository.keys())
        if not timestamps:
            return [], None

        # Iterate through timestamps in reverse order to find the most recent trajectory
        for timestamp in reversed(timestamps):
            state_data = data_repository.get(timestamp, {})
            flight_data = state_data.get(flight_key, {})
            rtg = flight_data.get("Flight", {}).get("routeTrajectoryGroup", {})
            predicted_points = rtg.get("predicted", {})
            if predicted_points:
                return predicted_points.get("element", []), timestamp

        return [], None

    def _filter_most_recent_predicted_trajectory(self, last_predicted_points, found_timestamp, current_timestamp):
        """
        Filters the most recent predicted trajectory points based on the found timestamp.
        If current_timestamp is provided, it filters the points to include only those after the current_timestamp.
        """

        # Ensure datetimes
        if isinstance(found_timestamp, str):
            found_timestamp = datetime.fromisoformat(found_timestamp)
        if isinstance(current_timestamp, str):
            current_timestamp = datetime.fromisoformat(current_timestamp) if current_timestamp else None

        # Round to seconds
        found_timestamp = found_timestamp.replace(microsecond=0)
        if current_timestamp:
            current_timestamp = current_timestamp.replace(microsecond=0)

        # No reference timestamp → return everything
        if not current_timestamp:
            return last_predicted_points

        # Calculate index directly
        delta_seconds = (current_timestamp - found_timestamp).total_seconds()
        start_idx = int(delta_seconds) - 1  # -1 because first point is +1s from found_timestamp

        # If the calculated index is out of range, prediction is outdated
        if start_idx < 0 or start_idx >= len(last_predicted_points):
            return []

        return last_predicted_points[start_idx:]   

    def _extract_repo_cd_data(self, data_repository, timestamp):
        """
        Extracts data from the data repository for conflict detection (CD).
        Gathers current position, flight level, predicted trajectory, cruise level,
        waypoints, speed, airport and within_sector information for each flight at the given timestamp.
        """
        state_data = data_repository.get(timestamp, {})
        cd_data = {timestamp: []}

        for flight_key, flight_data in state_data.items():
            # Ensure essential keys exist before accessing
            if not isinstance(flight_data, dict):
                return
            flight = flight_data.get("Flight")
            if not isinstance(flight, dict):
                return
            rtg = flight.get("routeTrajectoryGroup")
            if not isinstance(rtg, dict):
                return

            predicted_points = []
            # Conflict detection requires at least two flights to compare
            if len(state_data.keys()) >= 2:
                # Tražimo podskup prediktane trajektorije od trenutnog timestampa (zadanog kroz input) pa nadalje.
                # To je potrebno za sve avione, nekad će timestamp input biti početak trajetorije za neki avion, 
                # a za ostale ćemo raditi odrezivanje liste prediktane trajektorije.
                if "predicted" not in rtg or ("predicted" in rtg and not isinstance(rtg["predicted"], dict)):
                    last_predicted_points, found_timestamp = self._find_most_recent_predicted_trajectory(data_repository, flight_key)
                    if last_predicted_points:
                        predicted_points = self._filter_most_recent_predicted_trajectory(last_predicted_points, found_timestamp, timestamp)
                        # print("Original trajectory:", last_predicted_points)
                        # print("Filtered trajectory:", predicted_points)
                        if not predicted_points:
                            pass
                            # TODO: run our own trajectory prediction algorithm, since the most recent predicted trajectory is outdated
                    else:
                        pass
                        # TODO: run our own trajectory prediction algorithm, since there is no predicted trajectory available from Polaris
                else:
                    predicted_points = rtg["predicted"].get("element", [])
                    # print("Keeping whole trajectory: ", predicted_points)

            current_branch = rtg.get("current", {})
            if isinstance(current_branch, dict):
                current_point = current_branch.get("element", {}).get("point4D", {})
                current_pos = current_point.get("position", {}).get("pos")
                current_fl = current_point.get("level", {}).get("flightLevel")
                current_velocity = current_point.get("velocity", {})
            else:
                current_pos = current_fl = current_velocity = None

            cleared_branch = rtg.get("agreed", [])
            if not isinstance(cleared_branch, list):
                return
            cleared_waypoints = [
                entry.get("element", {}).get("point4D", {})
                for entry in cleared_branch
                if isinstance(entry, dict)
            ]
            if not cleared_waypoints:
                return

            cleared_route_info = next(
                (item for item in cleared_branch if isinstance(item, dict) and "routeInformation" in item),
                {}
            )
            cruise_level_cleared = (
                cleared_route_info.get("routeInformation", {})
                                .get("cruisingLevel", {})
                                .get("flightLevel")
            )
            if cruise_level_cleared is None:
                return

            desired_branch = rtg.get("desired", [])
            if not isinstance(desired_branch, list):
                return
            desired_route_info = next(
                (item for item in desired_branch if isinstance(item, dict) and "routeInformation" in item),
                {}
            )
            cruise_level_desired = (
                desired_route_info.get("routeInformation", {})
                                .get("cruisingLevel", {})
                                .get("flightLevel")
            )
            if cruise_level_desired is None:
                return

            exit_point = flight.get("enRoute", {}).get("exitPoint")
            departure_indicator = flight.get("departure", {}).get("departureAerodrome", {}).get("locationIndicator")

            arrival = flight.get("arrival", {})
            arrival_indicator = arrival.get("destinationAerodrome", {}).get("locationIndicator")

            alt_arrival = arrival.get("alternateAerodromeAlternate", {})
            alt_arrival_indicator = alt_arrival.get("locationIndicator") if isinstance(alt_arrival, dict) else None

            within_sector = flight.get("withinSector", None)

            cd_data[timestamp].append({
                "callsign": flight_key,
                "current_pos": current_pos,
                "current_fl": current_fl,
                "current_velocity": current_velocity,
                "predicted_points": predicted_points,
                "cruise_level_cleared": cruise_level_cleared,
                "cruise_level_desired": cruise_level_desired,
                "exit_point": exit_point,
                "cleared_waypoints": cleared_waypoints,
                "departure_aerodrome_indicator": departure_indicator,
                "arrival_aerodrome_indicator": arrival_indicator,
                "arrival_aerodrome_alternate_indicator": alt_arrival_indicator,
                "within_sector": within_sector
            })

        return cd_data

    def _find_conflicts(self, cd_data):
        """
        This method implements the logic to identify conflicts based on the data from repository.
        """
        pass
    
    def detect(self, data_repository, timestamp):
        """
        Detect conflicts in the data repository based on the timestamp.
        This method should handle data extraction and calling the main method for finding the conflicts.
        """

        cd_data = self._extract_repo_cd_data(data_repository, timestamp)
        # print("CD data extracted: ", cd_data)

        """conflicts = self._find_conflicts(cd_data)
        if conflicts:
            self.detections[timestamp] = conflicts"""
        return cd_data
    


class ConflictResolution:
    def __init__(self):
        self.resolutions = {} # key: timestamp, value: list of resolutions for conflicts

    def resolve(self, conflict_detections):
        """
        Resolve conflicts based on the detected conflicts.
        This method implement the logic to resolve conflicts.
        """
        pass