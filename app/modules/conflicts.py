import json
import re
from math import radians, sin, cos, atan2, sqrt
from datetime import datetime, timedelta
from typing import List, Dict, Any, Optional, Tuple
import math
import numpy as np
from dateutil import parser
from geopy.distance import geodesic

class ConflictPreprocessor:
    """
    Preprocess CD data (from data_repository).
    Input: cd_data produced by _extract_repo_cd_data (dict: timestamp -> list[flight_dicts])
    """

    TIMESTAMP_PATTERN = re.compile(r"Timestamp\('([^']+)'\)")

    def __init__(self, cd_data: Dict[str, List[Dict[str, Any]]]):
        # cd_data: { timestamp_str_or_dt: [ {callsign:..., predicted_points:..., cleared_waypoints:..., pred_times: ...}, ... ] }
        # Normalize keys to datetime objects for convenience
        self.cd_data = {}
        for k, v in cd_data.items():
            if isinstance(k, str):
                try:
                    key_dt = datetime.fromisoformat(k)
                except Exception:
                    key_dt = k  # keep as-is if cannot parse
            else:
                key_dt = k
            self.cd_data[key_dt] = v

    # ---------------------------
    # Utility functions
    # ---------------------------
    @staticmethod
    def _haversine_m(lat1, lon1, lat2, lon2) -> float:
        R = 6371000.0
        phi1, phi2 = radians(lat1), radians(lat2)
        dphi = radians(lat2 - lat1)
        dlambda = radians(lon2 - lon1)
        a = sin(dphi / 2)**2 + cos(phi1) * cos(phi2) * sin(dlambda / 2)**2
        c = 2 * atan2(sqrt(a), sqrt(1 - a))
        return R * c

    @classmethod
    def _parse_predicted_points(cls, pts_raw) -> Tuple[List[float], List[float], List[Any], List[Any]]:
        """
        pts_raw may be:
         - list of dicts: [{'latDeg':..., 'lonDeg':..., 'altitudeFt':..., 'trueHeading':...}, ...]
         - string (JSON) -> parse
         - None -> return empty lists
        returns (lats, lons, fls, headings)
        """
        if pts_raw is None:
            return [], [], [], []
        if isinstance(pts_raw, str):
            try:
                pts = json.loads(pts_raw)
            except Exception:
                return [], [], [], []
        elif isinstance(pts_raw, list):
            pts = pts_raw
        else:
            # unknown type
            return [], [], [], []

        lats, lons, fls, heads = [], [], [], []
        for p in pts:
            if not isinstance(p, dict):
                continue
            lats.append(p.get("latDeg"))
            lons.append(p.get("lonDeg"))
            # altitude might be named altitudeFt or flightLevel etc. keep both tries
            fls.append(p.get("altitudeFt") if "altitudeFt" in p else p.get("flightLevel"))
            heads.append(p.get("trueHeading"))
        return lats, lons, fls, heads

    @classmethod
    def _extract_cleared_coords(cls, cleared_waypoints_raw) -> List[Tuple[float, float]]:
        """
        cleared_waypoints_raw expected list of point4D dicts or JSON.
        Return list of (lat, lon) tuples.
        """
        coords = []
        if cleared_waypoints_raw is None:
            return coords
        if isinstance(cleared_waypoints_raw, str):
            try:
                cleared_waypoints = json.loads(cleared_waypoints_raw)
            except Exception:
                return coords
        else:
            cleared_waypoints = cleared_waypoints_raw
        if not isinstance(cleared_waypoints, list):
            return coords
        for entry in cleared_waypoints:
            if not isinstance(entry, dict):
                continue
            # entry expected to be point4D dict or wrapper with 'element'
            point4d = entry if "position" in entry else entry.get("element", {}).get("point4D", entry.get("point4D", {}))
            if not isinstance(point4d, dict):
                continue
            pos = point4d.get("position", {}) or {}
            pos = pos.get("pos", {}) if isinstance(pos, dict) else {}
            lat = pos.get("lat")
            lon = pos.get("lon")
            if lat is not None and lon is not None:
                coords.append((lat, lon))
        return coords

    @classmethod
    def _parse_pred_times(cls, times_field) -> List[datetime]:
        """
        times_field can be:
         - list of ISO timestamps (strings)
         - string containing "Timestamp('...')" occurrences
         - None
        Return list of datetimes (rounded to seconds).
        """
        if times_field is None:
            return []
        if isinstance(times_field, list):
            out = []
            for item in times_field:
                if not item:
                    continue
                try:
                    dt = datetime.fromisoformat(item)
                    out.append(dt.replace(microsecond=0))
                except Exception:
                    # try extracting from any string
                    if isinstance(item, str):
                        m = cls.TIMESTAMP_PATTERN.findall(item)
                        for mm in m:
                            try:
                                out.append(datetime.fromisoformat(mm).replace(microsecond=0))
                            except Exception:
                                pass
            return out
        if isinstance(times_field, str):
            matches = cls.TIMESTAMP_PATTERN.findall(times_field)
            out = []
            for m in matches:
                try:
                    out.append(datetime.fromisoformat(m).replace(microsecond=0))
                except Exception:
                    pass
            return out
        # fallback: try to stringify and extract
        try:
            text = str(times_field)
            matches = cls.TIMESTAMP_PATTERN.findall(text)
            out = []
            for m in matches:
                try:
                    out.append(datetime.fromisoformat(m).replace(microsecond=0))
                except Exception:
                    pass
            return out
        except Exception:
            return []

    @staticmethod
    def _initial_speed_kmh(lats: List[float], lons: List[float]) -> Optional[float]:
        if not lats or not lons or len(lats) < 2 or len(lons) < 2:
            return None
        if lats[0] is None or lons[0] is None or lats[1] is None or lons[1] is None:
            return None
        try:
            dist_m = ConflictPreprocessor._haversine_m(lats[0], lons[0], lats[1], lons[1])
            return dist_m * 3.6
        except Exception:
            return None

    # ---------------------------
    # Trajectory trimming (fast index calculation)
    # ---------------------------
    @staticmethod
    def _find_index_for_current(found_timestamp: datetime, current_timestamp: datetime) -> Optional[int]:
        """Return index in predicted list corresponding to current_timestamp. Rules:
           - first predicted point is found_timestamp + 1s (index 0 => found_ts + 1)
           - returns index or None if outside range
        """
        if found_timestamp is None or current_timestamp is None:
            return None
        # round to seconds
        f = found_timestamp.replace(microsecond=0)
        c = current_timestamp.replace(microsecond=0)
        delta = (c - f).total_seconds()
        idx = int(delta) - 1  # because index 0 corresponds to +1s
        return idx if idx >= 0 else None

    def _filter_predicted_trajectory(self, last_predicted_points, found_timestamp, current_timestamp):
        # TODO: refactor this since it is identical to ConflictDetection._filter_most_recent_predicted_trajectory(...)
        """
        last_predicted_points: list (raw)
        found_timestamp/current_timestamp: datetime or string
        returns trimmed list (points from matching index onwards) or [] if not match
        """
        if not last_predicted_points:
            return []
        if isinstance(found_timestamp, str):
            try:
                found_ts = datetime.fromisoformat(found_timestamp)
            except Exception:
                return last_predicted_points
        else:
            found_ts = found_timestamp
        if isinstance(current_timestamp, str):
            try:
                curr_ts = datetime.fromisoformat(current_timestamp)
            except Exception:
                curr_ts = None
        else:
            curr_ts = current_timestamp

        # if no current timestamp: keep all
        if curr_ts is None:
            return last_predicted_points

        idx = self._find_index_for_current(found_ts, curr_ts)
        if idx is None:
            return []
        if idx < 0:
            return []
        if idx >= len(last_predicted_points):
            return []
        return last_predicted_points[idx:]

    # ---------------------------
    # Preprocessing single timestamp
    # ---------------------------
    def preprocess_for_timestamp(self, timestamp_key) -> List[Dict[str, Any]]:
        """
        Preprocess flights for a given timestamp (key should be datetime or ISO string).
        Returns list of processed flight dicts ready for conflict detection.
        """
        if isinstance(timestamp_key, str):
            try:
                ts = datetime.fromisoformat(timestamp_key)
            except Exception:
                ts = timestamp_key
        else:
            ts = timestamp_key

        raw_list = self.cd_data.get(ts, [])
        if not raw_list:
            return []

        # Deduplicate by callsign (keep first), preserving order
        seen = set()
        deduped = []
        for rec in raw_list:
            callsign = rec.get("callsign")
            if callsign is None:
                continue
            if callsign in seen:
                continue
            seen.add(callsign)
            deduped.append(rec)

        processed = []
        for rec in deduped:
            callsign = rec.get("callsign")
            # predicted points may be in rec['predicted_points'] or rec['predicted_points'] as JSON string
            raw_pred = rec.get("predicted_points")
            # No predicted points, cannot process
            if not raw_pred:
                return
            
            # print("Processing flight:", callsign, "at timestamp:", ts)
            pred_lats, pred_lons, pred_fl, pred_heads = self._parse_predicted_points(raw_pred)

            # raw_pred_times = rec.get("predicted_times")
            # print("Raw predicted times:", len(raw_pred_times))
            # pred_times = self._parse_pred_times(raw_pred_times)
            # print("Len of parsed predicted times:", len(pred_times))
            pred_times = rec.get("predicted_times") # Already processed in _extract_repo_cd_data
            # print("Predicted times len:", len(pred_times))

            # If this flight contains a separate "found_timestamp" for that predicted trajectory,
            # it can be passed in rec or inferred. We'll try rec.get('found_timestamp') else ts.
            found_ts = rec.get("found_timestamp", ts)

            # Trim predicted trajectory to only points after current timestamp (ts)
            trimmed_predicted = self._filter_predicted_trajectory(raw_pred, found_ts, ts)
            # reparse trimmed predicted points lists to get lat/lon etc for trimmed list
            tr_lats, tr_lons, tr_fl, tr_heads = self._parse_predicted_points(trimmed_predicted)

            # cleared waypoints coords
            cleared_coords = self._extract_cleared_coords(rec.get("cleared_waypoints"))

            # initial speed
            speed_kmh = self._initial_speed_kmh(tr_lats if tr_lats else pred_lats, tr_lons if tr_lons else pred_lons)

            # map cruise fields to the names used in notebook if needed
            cod_flightlevel = rec.get("cruise_level_cleared") or rec.get("cod_flightlevel") or rec.get("cruise_level")
            ind_cruise_level = rec.get("cruise_level_desired") or rec.get("ind_cruise_level")

            processed.append({
                "timestamp": ts,
                "callsign": callsign,
                "current_pos": rec.get("current_pos"),
                "current_fl": rec.get("current_fl"),
                "current_velocity": rec.get("current_velocity"),
                "cod_flightlevel": cod_flightlevel,
                "ind_cruise_level": ind_cruise_level,
                "exit_point": rec.get("exit_point"),
                "cleared_waypoints": rec.get("cleared_waypoints"),
                "cleared_coords": cleared_coords,
                "departure_aerodrome_indicator": rec.get("departure_aerodrome_indicator"),
                "arrival_aerodrome_indicator": rec.get("arrival_aerodrome_indicator"),
                "arrival_aerodrome_alternate_indicator": rec.get("arrival_aerodrome_alternate_indicator"),
                "within_sector": rec.get("within_sector"),
                "predicted_points_raw": raw_pred,
                "predicted_points_trimmed": trimmed_predicted,
                "pred_lats": tr_lats if tr_lats else pred_lats,
                "pred_lons": tr_lons if tr_lons else pred_lons,
                "pred_fl": tr_fl if tr_fl else pred_fl,
                "pred_heads": tr_heads if tr_heads else pred_heads,
                "pred_times": pred_times,
                "initial_speed_kmh": speed_kmh
            })

        return processed

    # ---------------------------
    # Align / trim list-like fields (recortar_y_alinear_columnas)
    # ---------------------------
    @staticmethod
    def align_and_trim_rows(rows: List[Dict[str, Any]], list_keys: List[str], max_len: int = 1000) -> List[Dict[str, Any]]:
        """
        For each dict (row), ensure listed list_keys are trimmed to the same minimum length (but <= max_len).
        Rows where lists are missing are left as-is.
        Returns new list of dicts (modified in place).
        """
        for row in rows:
            # find min length among present lists
            lens = [len(row.get(k, [])) for k in list_keys if isinstance(row.get(k, []), list)]
            if not lens:
                continue
            min_len = min(min(lens), max_len)
            for k in list_keys:
                if isinstance(row.get(k), list):
                    row[k] = row[k][:min_len]
        return rows

    # ---------------------------
    # Split into contiguous segments <= interval minutes (no roman numerals)
    # ---------------------------
    @staticmethod
    def split_into_segments(rows: List[Dict[str, Any]], interval_minutes: int = 5) -> List[List[Dict[str, Any]]]:
        if not rows:
            return []
        # assume rows have 'timestamp' as datetime
        rows_sorted = sorted(rows, key=lambda r: r.get("timestamp"))
        segments = []
        start = 0
        n = len(rows_sorted)
        while start < n:
            start_ts = rows_sorted[start]["timestamp"]
            end = start
            while end + 1 < n and (rows_sorted[end + 1]["timestamp"] - start_ts) <= timedelta(minutes=interval_minutes):
                end += 1
            segments.append(rows_sorted[start:end + 1])
            start = end + 1
        return segments

    # ---------------------------
    # Top-level API to get sample for CD algorithm
    # ---------------------------
    def get_preprocessed_sample(self, timestamp_key, dedupe: bool = True, list_keys: Optional[List[str]] = None,
                                align_max_len: int = 1000) -> List[Dict[str, Any]]:
        """
        Returns a list of processed flight dicts for the given timestamp,
        deduplicated and trimmed/aligned ready for conflict detection.
        """
        rows = self.preprocess_for_timestamp(timestamp_key)
        # Optionally dedupe already done in preprocess_for_timestamp
        # Align/trim list fields
        if rows and list_keys:
            rows = self.align_and_trim_rows(rows, list_keys, max_len=align_max_len)
            return rows
        return



class ConflictDetection:
    def __init__(self):
        """
        detections: ...
        """
        self.detections = {}

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

    
    def _generate_times_from_trajectory(self, predicted_points, timestamp):
        """
        Generates a list of timestamps for each predicted point based on the initial timestamp.
        Assumes predicted_points is a list of dicts.
        
        Args:
            predicted_points (list): List of dicts, one per predicted point.
            timestamp (str or datetime): Initial timestamp (ISO string or datetime).
        
        Returns:
            list[datetime]: List of timestamps, one for each predicted point.
        """
        if not predicted_points:
            return []

        # Ensure timestamp is a datetime object
        if isinstance(timestamp, str):
            timestamp = datetime.fromisoformat(timestamp)
        
        timestamp = timestamp.replace(microsecond=0)  # Round to seconds

        # Generate timestamps with 1-second intervals
        trajectory_len = len(predicted_points)
        times = [timestamp + timedelta(seconds=i) for i in range(1, trajectory_len + 1)]

        return times

    def _extract_repo_cd_data(self, data_repository, state_data, timestamp):
        """
        Extracts data from the data repository for conflict detection (CD).
        Gathers current position, flight level, predicted trajectory, cruise level,
        waypoints, speed, airport and within_sector information for each flight at the given timestamp.
        """

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
            
            # Tražimo podskup prediktane trajektorije od trenutnog timestampa (zadanog kroz input) pa nadalje.
            # To je potrebno za sve avione, nekad će timestamp input biti početak trajetorije za neki avion, 
            # a za ostale ćemo raditi odrezivanje liste prediktane trajektorije.
            if "predicted" not in rtg or ("predicted" in rtg and not isinstance(rtg["predicted"], dict)):
                last_predicted_points, found_timestamp = self._find_most_recent_predicted_trajectory(data_repository, flight_key)
                if last_predicted_points:
                    predicted_points = self._filter_most_recent_predicted_trajectory(last_predicted_points, found_timestamp, timestamp)
                    # print("FLight key: ", flight_key, "timestamp: ", timestamp)
                    # print("Original trajectory len: ", len(last_predicted_points))
                    # print("Filtered trajectory len: ", len(predicted_points))
                    if not predicted_points:
                        pass
                        # TODO: run our own trajectory prediction algorithm, since the most recent predicted trajectory is outdated
                else:
                    pass
                    # TODO: run our own trajectory prediction algorithm, since there is no predicted trajectory available from Polaris
            else:
                predicted_points = rtg["predicted"].get("element", [])
                # print("Keeping whole trajectory: ", predicted_points)

            predicted_times = self._generate_times_from_trajectory(predicted_points, timestamp)

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
                cleared_route_info.get("routeInformation", {}).get("cruisingLevel", {}).get("flightLevel")
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
                "current_velocity": current_velocity,   # TODO: this should be speed
                "predicted_points": predicted_points,
                "predicted_times": predicted_times,
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
    
    # ---------------------------
    # Utility functions
    # ---------------------------

    def convert_to_timestamp(self, time_obj):
        """
        Converts any time object (str, datetime, int, float) to naive datetime (no timezone).
        """
        if time_obj is None:
            return None

        try:
            # Handle NaN safely (for floats)
            if isinstance(time_obj, float) and math.isnan(time_obj):
                return None

            if isinstance(time_obj, str):
                time_obj = parser.parse(time_obj)
            elif isinstance(time_obj, (int, float)):  # treat as UNIX timestamp
                time_obj = datetime.fromtimestamp(time_obj)
            elif not isinstance(time_obj, datetime):
                # Unknown type
                return None

            # Ensure naive datetime (strip timezone if present)
            return time_obj.replace(tzinfo=None) if time_obj.tzinfo else time_obj

        except Exception as e:
            print(f"Error in convert_to_timestamp: {e}")
            return None
    
    def _hms_to_seconds(self, hms: str) -> int:
        # Parse as a time object
        t = datetime.strptime(hms, "%H:%M:%S").time()
        td = timedelta(hours=t.hour, minutes=t.minute, seconds=t.second)
        return int(td.total_seconds())


    def calculate_intermediate_point(self, lat1, lon1, lat2, lon2):
        """
        Calculates intermediate point between two geographic coordinates using spherical navigation formulas.
        Returns latitude and longitude in degrees.
        """
        lat1, lon1 = np.radians(lat1), np.radians(lon1)
        lat2, lon2 = np.radians(lat2), np.radians(lon2)
        dlon = lon2 - lon1

        Bx = np.cos(lat2) * np.cos(dlon)
        By = np.cos(lat2) * np.sin(dlon)

        lat_mid = np.arctan2(np.sin(lat1) + np.sin(lat2),
                            np.sqrt((np.cos(lat1) + Bx)**2 + By**2))
        lon_mid = lon1 + np.arctan2(By, np.cos(lat1) + Bx)

        return np.degrees(lat_mid), np.degrees(lon_mid)

    """def format_time(self, x):
        try:
            import pandas as pd
            if pd.isna(x):
                return None
            ts = self.convert_to_timestamp(x)
            if ts is not None:
                return ts.strftime('%H:%M:%S')
            return None
        except Exception:
            return None"""
    
    def format_time(self, x):
        try:
            if x is None:
                return None
            ts = self.convert_to_timestamp(x)
            if ts is not None:
                return ts.strftime('%H:%M:%S')
            return None
        except Exception:
            return None


    def calculate_cumulative_distances(self, latitudes, longitudes, index_min_dis):
        """
        Calculate cumulative geodesic distances in nautical miles between consecutive points up to index_min_dis.
        """
        if not isinstance(index_min_dis, int) or index_min_dis >= len(latitudes):
            return None
        
        distances = []
        for i in range(index_min_dis):
            pt1 = (latitudes[i], longitudes[i])
            pt2 = (latitudes[i + 1], longitudes[i + 1])
            dist_nm = geodesic(pt1, pt2).nautical
            distances.append(dist_nm)
        return sum(distances)
    
    # ---------------------------
    # CD algorithm functions
    # ---------------------------
    
    def generate_aircraft_pairs_bidirectional(self, sample):
        """
        Generates all pairs (A,B) and (B,A) excluding (A,A) from sample.
        Each pair is a combined dict with suffixed keys.
        Returns list of dicts representing pairs.
        """
        callsigns = [ac['callsign'] for ac in sample]
        pairs = [(a1, a2) for a1 in callsigns for a2 in callsigns if a1 != a2]
        results = []
        
        # Create a dict from callsign to aircraft data for fast lookup
        ac_map = {ac['callsign']: ac for ac in sample}

        for ac1, ac2 in pairs:
            row_ac1 = {f"{k}_1": v for k, v in ac_map[ac1].items()}
            # print("ac1 data:", len(row_ac1["pred_lats_1"]), len(row_ac1["pred_lons_1"]), len(row_ac1["pred_fl_1"]), len(row_ac1["pred_times_1"]))
            row_ac2 = {f"{k}_2": v for k, v in ac_map[ac2].items()}
            # print("ac2 data:", len(row_ac2["pred_lats_2"]), len(row_ac2["pred_lons_2"]), len(row_ac2["pred_fl_2"]), len(row_ac2["pred_times_2"]))
            combined = {**row_ac1, **row_ac2}
            combined['Aircraft Pair'] = (ac1, ac2)
            results.append(combined)
        
        # for pair in results:
            # print("Aircraft pair generated:", pair['Aircraft Pair'])
        return results

    def separation_with_time_matching(self, pairs):
        """
        Adds separation info (horizontal, vertical, combined) for each pair dict in pairs list.
        Only times that match exactly (rounded to seconds) are considered.
        """
        def haversine(pos1, pos2):
            # Haversine formula returning nautical miles
            R = 3440.065  # Earth radius in nautical miles
            lat1, lon1 = np.radians(pos1[0]), np.radians(pos1[1])
            lat2, lon2 = np.radians(pos2[0]), np.radians(pos2[1])
            dlat = lat2 - lat1
            dlon = lon2 - lon1
            a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
            c = 2 * np.arcsin(np.sqrt(a))
            return R * c
        
        for row in pairs:
            horiz_sep = []
            vert_sep = []
            combined_sep = []

            lat1 = row['pred_lats_1']
            # print("Length of lat1:", len(lat1))
            lon1 = row['pred_lons_1']
            # print("Length of lon1:", len(lon1))
            fl1 = row['pred_fl_1']
            # print("Length of fl1:", len(fl1))
            # print("Length of raw times1:", len(row['pred_times_1']))
            # times1 = [self.convert_to_timestamp(t).replace(microsecond=0) for t in row['pred_times_1']]
            times1 = row['pred_times_1']  # Already processed in _extract_repo_cd_data
            # print("Length of times1:", len(times1))

            lat2 = row['pred_lats_2']
            # print("Length of lat2:", len(lat2))
            lon2 = row['pred_lons_2']
            # print("Length of lon2:", len(lon2))
            fl2 = row['pred_fl_2']
            # print("Length of fl2:", len(row['pred_fl_2']))
            # print("Length of raw times2:", len(row['pred_times_2']))
            # times2 = [self.convert_to_timestamp(t).replace(microsecond=0) for t in row['pred_times_2']]
            times2 = row['pred_times_2']  # Already processed in _extract_repo_cd_data
            # print("Length of times2:", len(times2))

            for i, t1 in enumerate(times1):
                try:
                    j = times2.index(t1)
                except ValueError:
                    print("Time", t1, "not found in times2 for pair:", row['Aircraft Pair'])
                    continue
                d_horiz = haversine((lat1[i], lon1[i]), (lat2[j], lon2[j]))
                # print("Horizontal distance at time", t1, ":", d_horiz)
                d_vert = fl1[i] - fl2[j]
                # print("Vertical distance at time", t1, ":", d_vert)
                horiz_sep.append((d_horiz, t1))
                vert_sep.append((d_vert, t1))
                combined_sep.append((d_horiz, d_vert, t1))

            row['horizontal_separation'] = horiz_sep
            row['vertical_separation'] = vert_sep
            row['separation'] = combined_sep

            # print("Separation data added to pairs:", row["horizontal_separation"], row["vertical_separation"], row["separation"])
        return pairs

    def process_conflicts_and_SI(self, pairs):
        """
        For each pair dict, sets 'SI' and 'Conflict' flags based on separation data.
        SI if horizontal < 10 nm and vertical < 1000 ft.
        Conflict if horizontal < 5 nm and vertical < 1000 ft.
        """
        for row in pairs:
            si_condition = False
            conflict_condition = False

            for h_sep, v_sep, _ in row.get('separation', []):
                if h_sep < 10 and abs(v_sep) < 1000:
                    si_condition = True
                    if h_sep < 5:
                        conflict_condition = True
                        break
                elif h_sep < 5 and abs(v_sep) < 1000:
                    conflict_condition = True
                    si_condition = True
                    break
            
            row['SI'] = 1 if si_condition else 0
            row['Conflict'] = 1 if conflict_condition else 0
        
        # print("Processed conflicts and SI flags:", pairs)
        return pairs

    def process_separations_MinSep_and_tcpa_ONLY_SI(self, pairs):
        """
        For pairs with SI==1, calculates conditioned separations, combined separations,
        minimum combined separation value, TCPA, and min separation distances.
        Adds these to each pair dict.
        """
        for row in pairs:
            if row.get('SI', 0) == 1:
                separations = row.get('separation', [])
                conditioned_sep = [
                    (h, v, t) for h, v, t in separations if h < 10 and abs(v) < 1000
                ]
                if not conditioned_sep:
                    conditioned_sep = None
                
                if conditioned_sep:
                    combined_values = [h + (abs(v)/1000)*5 for h,v,_ in conditioned_sep]
                else:
                    combined_values = None

                if combined_values:
                    min_value = min(combined_values)
                    min_index = combined_values.index(min_value)
                    time_min = conditioned_sep[min_index][2]
                else:
                    min_value = None
                    min_index = None
                    time_min = None

                if conditioned_sep and min_index is not None:
                    h, v, t = conditioned_sep[min_index]
                    h_m = h * 1852  # nm to meters
                    v_m = v * 0.3048  # ft to meters
                    modulo_distance = math.sqrt(h_m**2 + v_m**2)
                    min_dis = (h, v, t, modulo_distance)
                else:
                    min_dis = None
            else:
                conditioned_sep = None
                combined_values = None
                min_value = None
                time_min = None
                min_dis = None
            
            row['conditioned_separation'] = conditioned_sep
            row['combined_separation'] = combined_values
            row['min_combined_sep'] = min_value
            row['tcpa'] = time_min
            row['min_sep'] = min_dis
        
        return pairs

    def process_cpa_and_times_ONLY_SI(self, pairs):
        """
        For pairs with SI==1, computes CPA positions for both aircraft and intermediate point,
        time to CPA from start for each aircraft, and indexes of CPA points.
        """
        for row in pairs:
            if row.get('SI', 0) != 1:
                row['position_cpa_1st_aircraft'] = None
                row['position_cpa_2nd_aircraft'] = None
                row['position_cpa_intermediate'] = None
                row['tcpa_hms'] = None
                row['time_to_cpa_hms_1'] = None
                row['time_to_cpa_hms_2'] = None
                row['index_cpa_1st_aircraft'] = None
                row['index_cpa_2nd_aircraft'] = None
                continue
            
            try:
                # tcpa = self.convert_to_timestamp(row['tcpa'])
                tcpa = row["tcpa"]
                if tcpa is None:
                    raise ValueError("Invalid TCPA")
                # times_1 = [self.convert_to_timestamp(t) for t in row['pred_times_1']]
                # times_2 = [self.convert_to_timestamp(t) for t in row['pred_times_2']]
                times_1 = row['pred_times_1']
                times_2 = row['pred_times_2']

                index_2 = min(range(len(times_2)), key=lambda i: abs(times_2[i] - tcpa))
                try:
                    index_1 = min(range(len(times_1)), key=lambda i: abs(times_1[i] - tcpa))
                except Exception:
                    index_1 = index_2

                pos_1 = (row['pred_lats_1'][index_1], row['pred_lons_1'][index_1], row['pred_fl_1'][index_1], times_1[index_1])
                pos_2 = (row['pred_lats_2'][index_2], row['pred_lons_2'][index_2], row['pred_fl_2'][index_2], times_2[index_2])

                lat_i, lon_i = self.calculate_intermediate_point(pos_1[0], pos_1[1], pos_2[0], pos_2[1])
                h_i = (pos_1[2] + pos_2[2]) / 2
                pos_intermediate = (lat_i, lon_i, h_i, tcpa)

                time_start_1 = times_1[0]
                delta_1 = tcpa - time_start_1
                ttcpa_1 = f"{int(delta_1.total_seconds()//3600):02d}:{int((delta_1.total_seconds()%3600)//60):02d}:{int(delta_1.total_seconds()%60):02d}"

                time_start_2 = times_2[0]
                delta_2 = tcpa - time_start_2
                ttcpa_2 = f"{int(delta_2.total_seconds()//3600):02d}:{int((delta_2.total_seconds()%3600)//60):02d}:{int(delta_2.total_seconds()%60):02d}"

                row['position_cpa_1st_aircraft'] = pos_1
                row['position_cpa_2nd_aircraft'] = pos_2
                row['position_cpa_intermediate'] = pos_intermediate
                row['tcpa_hms'] = self.format_time(tcpa)
                row['time_to_cpa_hms_1'] = ttcpa_1
                row['time_to_cpa_hms_2'] = ttcpa_2
                row['index_cpa_1st_aircraft'] = index_1
                row['index_cpa_2nd_aircraft'] = index_2
            except Exception as e:
                print(f"Error processing CPA and times: {e}")
                row['position_cpa_1st_aircraft'] = None
                row['position_cpa_2nd_aircraft'] = None
                row['position_cpa_intermediate'] = None
                row['tcpa_hms'] = None
                row['time_to_cpa_hms_1'] = None
                row['time_to_cpa_hms_2'] = None
                row['index_cpa_1st_aircraft'] = None
                row['index_cpa_2nd_aircraft'] = None

        return pairs

    def calculate_distance_to_cpa_ONLY_SI(self, pairs):
        """
        For pairs with SI==1, calculates cumulative distances to CPA for each aircraft.
        """
        for row in pairs:
            if row.get('SI', 0) != 1:
                row['distance_to_cpa_1st_aircraft'] = None
                row['distance_to_cpa_2nd_aircraft'] = None
                continue
            try:
                lats_1 = row['pred_lats_1']
                lons_1 = row['pred_lons_1']
                idx_1 = row['index_cpa_1st_aircraft']

                lats_2 = row['pred_lats_2']
                lons_2 = row['pred_lons_2']
                idx_2 = row['index_cpa_2nd_aircraft']

                dist1 = self.calculate_cumulative_distances(lats_1, lons_1, idx_1)
                dist2 = self.calculate_cumulative_distances(lats_2, lons_2, idx_2)

                row['distance_to_cpa_1st_aircraft'] = dist1
                row['distance_to_cpa_2nd_aircraft'] = dist2
            except Exception as e:
                print(f"Error calculating distance to CPA: {e}")
                row['distance_to_cpa_1st_aircraft'] = None
                row['distance_to_cpa_2nd_aircraft'] = None
        
        return pairs

    def process_SI_pairs(self, pairs):
        """
        Filters pairs where SI == 1 and removes duplicates by unordered pairs.
        """
        # Filter pairs with SI==1
        si_pairs = [row for row in pairs if row.get('SI', 0) == 1]

        # Remove duplicates based on sorted aircraft pair callsigns
        seen = set()
        unique_si_pairs = []
        for row in si_pairs:
            pair_key = tuple(sorted(row['Aircraft Pair']))
            if pair_key not in seen:
                seen.add(pair_key)
                unique_si_pairs.append(row)
        
        # Convert flight level fields to int if they exist
        for row in unique_si_pairs:
            for key in ['cod_flightlevel_1', 'ind_cruise_level_1', 'cod_flightlevel_2', 'ind_cruise_level_2']:
                if key in row and row[key] is not None:
                    row[key] = int(row[key]["flightLevel"])
        
        return unique_si_pairs

    def add_initial_status(self, pairs):
        """
        Adds 'Status_initial_1' and 'Status_initial_2' indicating climb/descend/cruise based on first two FLs in pred_fl lists.
        """
        def get_status(fl_list):
            if not isinstance(fl_list, (list, tuple)) or len(fl_list) < 2:
                return None
            if fl_list[0] > fl_list[1]:
                return 'Descend'
            elif fl_list[0] < fl_list[1]:
                return 'Climb'
            else:
                return 'Cruise'

        for row in pairs:
            row['Status_initial_1'] = get_status(row.get('pred_fl_1', []))
            row['Status_initial_2'] = get_status(row.get('pred_fl_2', []))
        return pairs

    def add_status_ONLY_SI(self, pairs):
        """
        Adds 'Status_CPA_1' and 'Status_CPA_2' based on altitude at CPA compared to cod_flightlevel * 100.
        """
        def get_altitude_safe(pos_tuple):
            if pos_tuple is None or len(pos_tuple) < 3:
                return float('nan')
            return pos_tuple[2]
        
        for row in pairs:
            if row.get('SI', 0) != 1:
                row['Status_CPA_1'] = ''
                row['Status_CPA_2'] = ''
                continue
            a1 = get_altitude_safe(row.get('position_cpa_1st_aircraft'))
            c1 = row.get('cod_flightlevel_1', None)
            a2 = get_altitude_safe(row.get('position_cpa_2nd_aircraft'))
            c2 = row.get('cod_flightlevel_2', None)

            if a1 is None or c1 is None:
                row['Status_CPA_1'] = 'Unknown'
            else:
                c1_m = c1 * 100
                if a1 < c1_m:
                    row['Status_CPA_1'] = 'Descend'
                elif a1 > c1_m:
                    row['Status_CPA_1'] = 'Climb'
                else:
                    row['Status_CPA_1'] = 'Cruise'
            
            if a2 is None or c2 is None:
                row['Status_CPA_2'] = 'Unknown'
            else:
                c2_m = c2 * 100
                if a2 < c2_m:
                    row['Status_CPA_2'] = 'Descend'
                elif a2 > c2_m:
                    row['Status_CPA_2'] = 'Climb'
                else:
                    row['Status_CPA_2'] = 'Cruise'

        return pairs

    def process_conflict_pairs(self, pairs):
        """
        Filter only conflict cases (Conflict == 1), remove duplicates by unordered callsign pair,
        and ensure numeric flight levels are integers.
        
        Args:
            pairs (list of dict): Each dict has keys like 'Conflict', 'callsign_1', 'callsign_2', 
                                'cod_flightlevel_1', etc.
        
        Returns:
            list of dict: Filtered conflict pairs with numeric flight level fields.
        """
        seen_pairs = set()
        unique_conflicts = []

        for row in pairs:
            if row.get('Conflict') != 1:
                continue  # skip if not a conflict
            
            # Create an unordered key for the aircraft pair
            pair_key = tuple(sorted([row.get('callsign_1'), row.get('callsign_2')]))
            if pair_key in seen_pairs:
                continue  # skip duplicates
            
            seen_pairs.add(pair_key)

            # Convert numeric fields to int (if possible)
            for col in ['cod_flightlevel_1', 'ind_cruise_level_1',
                        'cod_flightlevel_2', 'ind_cruise_level_2']:
                if col in row and row[col] is not None:
                    try:
                        row[col] = int(row[col])
                    except (ValueError, TypeError):
                        pass  # keep original if conversion fails

            # Filter only needed conflict information for the graph
            row_filtered = {
                "Aircraft Pair": row["Aircraft Pair"], 
                "Conflict": row["Conflict"],
                "position_cpa_1st_aircraft": row["position_cpa_1st_aircraft"],
                "position_cpa_2nd_aircraft": row["position_cpa_2nd_aircraft"],
                "position_cpa_intermediate": row["position_cpa_intermediate"],
                "time_to_cpa_sec_1": self._hms_to_seconds(row["time_to_cpa_hms_1"]),
                "time_to_cpa_sec_2": self._hms_to_seconds(row["time_to_cpa_hms_2"]),
                "distance_to_cpa_1st_aircraft": row["distance_to_cpa_1st_aircraft"],
                "distance_to_cpa_2nd_aircraft": row["distance_to_cpa_2nd_aircraft"],
            }
            
            unique_conflicts.append(row_filtered)

        return unique_conflicts

    
    def detect(self, data_repository, timestamp):
        """
        Detect conflicts in the data repository based on the timestamp.
        This method should handle data extraction and calling the main method for finding the conflicts.
        """

        state_data = data_repository.get(timestamp, {})

        # Conflict detection only if there are at least 2 aircraft in the state data
        if len(state_data.keys()) >= 2:
            # Step 1: Extract raw data for the timestamp (list of dicts, each dict = one aircraft)
            cd_data = self._extract_repo_cd_data(data_repository, state_data, timestamp)

            # Step 2: Preprocess data to get the sample in expected format (list of dicts)
            prep = ConflictPreprocessor(cd_data)
            sample = prep.get_preprocessed_sample(
                timestamp, list_keys=['pred_lats','pred_lons','pred_fl','pred_times']
            )

            if not sample:
                print(f"Skipping conflict detection, data not gathered yet for timestamp {timestamp}.")
                return

            # Step 3: Generate bidirectional aircraft pairs with suffixed keys (_1, _2)
            pairs = self.generate_aircraft_pairs_bidirectional(sample)

            # Step 4: Calculate separations where predicted times match
            pairs = self.separation_with_time_matching(pairs)

            # Step 5: Determine SI and Conflict flags based on separations
            pairs = self.process_conflicts_and_SI(pairs)

            # Step 6: Calculate combined separations, MinSep, TCPA (only where SI==1)
            pairs = self.process_separations_MinSep_and_tcpa_ONLY_SI(pairs)

            # Step 7: Calculate CPA positions, indexes, and times (only where SI==1)
            pairs = self.process_cpa_and_times_ONLY_SI(pairs)

            # Step 8: Calculate distance to CPA for each aircraft (only where SI==1)
            pairs = self.calculate_distance_to_cpa_ONLY_SI(pairs)

            # Step 9: Extract unique SI pairs (filter duplicates by unordered pair keys)
            si_pairs = self.process_SI_pairs(pairs)

            # Step 10: Add initial climb/cruise/descend status based on predicted FL trajectory
            si_pairs = self.add_initial_status(si_pairs)

            # Step 11: Add climb/cruise/descend status at CPA point (only SI pairs)
            si_pairs = self.add_status_ONLY_SI(si_pairs)

            # Step 12: Extract unique conflict pairs (filter duplicates by unordered pair keys, Conflict==1)
            conflicts = self.process_conflict_pairs(pairs)

            # Store and return results (store in detections dict per timestamp)
            self.detections[timestamp] = {
                'conflicts': conflicts
            }
            
            print("Conflicts: ", conflicts)
            return self.detections[timestamp]
        
        print(f"No conflicts detected at {timestamp}.")
        return
    


class ConflictResolution:
    def __init__(self):
        self.resolutions = {} # key: timestamp, value: list of resolutions for conflicts

    def resolve(self, conflict_detections):
        """
        Resolve conflicts based on the detected conflicts.
        This method implement the logic to resolve conflicts.
        """
        pass