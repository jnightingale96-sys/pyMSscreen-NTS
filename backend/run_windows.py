"""Native (no Docker) entry point used by START_pyMSscreen-NTS.bat: serves API + built frontend on
http://127.0.0.1:<PORT> with waitress and opens the browser. Bound to localhost only."""
import os
import sys
import threading
import webbrowser

if __name__ == "__main__":
    # the multiprocessing children spawned by the extraction route re-import this file; keep work under the guard
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from waitress import serve
    from web_app import app, BASE_WORKING_DIR, FRONTEND_DIST

    port = int(os.getenv("PORT", "8080"))
    url = f"http://127.0.0.1:{port}/"
    os.makedirs(BASE_WORKING_DIR, exist_ok=True)
    print(f"pyMSscreen-NTS  ->  {url}")
    print(f"  default working directory: {BASE_WORKING_DIR}")
    print(f"  frontend: {FRONTEND_DIST}")
    print("  Close this window to stop the application.")
    if os.getenv("NO_BROWSER") != "1":
        threading.Timer(1.5, lambda: webbrowser.open(url)).start()
    serve(app, host="127.0.0.1", port=port, threads=8, channel_timeout=3600)
