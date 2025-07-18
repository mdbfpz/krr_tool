import httpx
from rdfox_db.core import RDFoxDB
import os


class RDFoxQuery(RDFoxDB):
    def __init__(self, connection_id):
        """ Inherit selected properties from RDFoxDB and add specific ones for Query. """

        super().__init__()
        self.connection_id = connection_id
        
        self.headers_output_turtle = self.headers.copy()
        self.headers_output_turtle.update({"Accept": "text/turtle; charset=UTF-8"})
        self.headers_output_csv = self.headers.copy()
        self.headers_output_csv.update({"Accept": "text/csv; charset=UTF-8"})

        self.import_datalog_rules()
        self.verify_imported_rules()

    async def import_turtle_data(self, ttl_data: str):
        """Import Turtle data into the RDFox data store asynchronously, including prefixes."""
        url = (
            f"{self.endpoint}/datastores/{self.data_store}/content"
            f"?connection={self.connection_id}&operation=add-content-update-prefixes"
        )

        headers = {**self.headers}
        # Either let RDFox auto-detect or explicitly set Turtle type
        if "Content-Type" not in headers:
            headers["Content-Type"] = "text/turtle"

        async with httpx.AsyncClient() as client:
            try:
                print("Importing here:", url)
                resp = await client.patch(
                    url,
                    data=ttl_data,
                    headers=headers
                )
                if resp.status_code in self.valid_responses:
                    print("Import succeeded (prefixes included).")
                else:
                    print("Import failed with status:", resp.status_code)
                    print("Response:", resp.text)

            except httpx.HTTPStatusError as e:
                print(f"HTTP status error: {e.response.status_code} — {e.response.text}")
            except httpx.RequestError as e:
                print("Request error:", e)
            except Exception as e:
                print("Unexpected error during import:", e)

    def import_datalog_rules(self):
        """Synchronously import each Datalog (.dlog) rule file into RDFox."""
        script_dir = os.path.dirname(os.path.abspath(__file__))
        rules_folder = os.path.abspath(os.path.join(script_dir, "../src/rules"))

        dlog_files = [
            os.path.join(rules_folder, f)
            for f in os.listdir(rules_folder)
            if f.endswith(".dlog")
        ]

        if not dlog_files:
            print("No .dlog files found in", rules_folder)
            return

        url = (
            f"{self.endpoint}/datastores/{self.data_store}/content"
            f"?connection={self.connection_id}"
        )

        headers = {
            **self.headers,
            "Content-Type": "application/x.datalog",  # Correct MIME type for Datalog rules
        }

        print(f"Importing Datalog rules from: {rules_folder} ...")

        for filepath in dlog_files:
            filename = os.path.basename(filepath)
            try:
                with open(filepath, mode="r", encoding="utf-8") as f:
                    rules_content = f.read()

                with httpx.Client() as client:
                    resp = client.post(url, data=rules_content, headers=headers)
                    if resp.status_code in self.valid_responses:
                        print(f"[OK] Imported {filename}")
                    else:
                        print(f"[ERROR] Failed to import {filename}")
                        print("Status:", resp.status_code)
                        print("Response:", resp.text)

            except Exception as e:
                print(f"[EXCEPTION] Importing {filename} failed: {e}")

    def verify_imported_rules(self):
        """
        Check that Datalog rules have been imported into the datastore by
        retrieving its content (only rules).
        """
        url = (
            f"{self.endpoint}/datastores/{self.data_store}/content"
            f"?rule-domain=user&connection={self.connection_id}"
        )
        headers = {
            **self.headers,
            "Accept": "application/x.datalog"
        }

        print("Verifying if rules are present in the data store...")

        try:
            with httpx.Client() as client:
                resp = client.get(url, headers=headers)

                if resp.status_code in self.valid_responses:
                    text = resp.text.strip()
                    if text:
                        print("[OK] Retrieved rules.")
                        return text
                    else:
                        print("[WARN] No rules found in the 'user' domain.")
                        return ""
                else:
                    print("[ERROR] Could not retrieve datastore content")
                    print("Status:", resp.status_code)
                    print("Response:", resp.text)
                    return None

        except Exception as e:
            print(f"[EXCEPTION] Error while verifying rules: {e}")
            return None

    
    async def select_all(self):
        """Select all facts from the RDFox data store."""

        async with httpx.AsyncClient() as client:
            try:
                # Change query parameters based on select type
                import_resp = await client.get(
                    #f"{self.endpoint}/datastores/{self.data_store}/content?connection={self.connection_id}",
                    #f"{self.endpoint}/datastores/{self.data_store}/content?default",
                    f"{self.endpoint}/datastores/{self.data_store}/content",
                    headers=self.headers_output_turtle
                )
                # print("Fetched all from the datastore:\n", import_resp.content)
                if import_resp.status_code in self.valid_responses:
                    print("Query 'select_all' succeeded.")
                else:
                    print("Query 'select_all' failed.")
                return

            except httpx.HTTPStatusError as e:
                print(f"HTTP error occurred: {e.response.status_code} - {e.response.text}")
                return
            except httpx.RequestError as e:
                print(f"Request error occurred: {e}")
                return
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                return
    