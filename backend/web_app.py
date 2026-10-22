import os

from flask import Flask, jsonify, send_from_directory
from flask_cors import CORS

from app.config import BASE_WORKING_DIR, PORT, SECRET_KEY
from app.routes.state_routes import state_bp
from app.routes.table_routes import table_bp
from app.routes.extraction_routes import extraction_bp
from app.routes.data_routes import data_bp
from app.routes.pdf_export import register_pdf_export
from app.routes.nts_routes import nts_bp

# Built React app (frontend/dist). In Docker nginx serves it instead; natively on Windows Flask does.
FRONTEND_DIST = os.path.abspath(os.getenv(
    "FRONTEND_DIST", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "frontend", "dist")))
SPA_ROUTES = ["/", "/extraction", "/plotting", "/nts", "/semiquant", "/product-ions"]


app = Flask(__name__)
app.secret_key = SECRET_KEY

CORS(app, supports_credentials=True, origins="*")

register_pdf_export(app)

app.register_blueprint(state_bp)
app.register_blueprint(table_bp)
app.register_blueprint(extraction_bp)
app.register_blueprint(data_bp)
app.register_blueprint(nts_bp)


def _spa_index():
    if not os.path.exists(os.path.join(FRONTEND_DIST, "index.html")):
        return jsonify({"status": "error", "message": f"Frontend not built ({FRONTEND_DIST}). Run: npm run build in frontend/"}), 404
    return send_from_directory(FRONTEND_DIST, "index.html")


for _route in SPA_ROUTES:
    app.add_url_rule(_route, f"spa_{_route.strip('/') or 'root'}", _spa_index)


@app.route("/assets/<path:name>")
def _spa_assets(name):
    return send_from_directory(os.path.join(FRONTEND_DIST, "assets"), name)


@app.route("/images/<path:name>")
def _spa_images(name):
    return send_from_directory(os.path.join(FRONTEND_DIST, "images"), name)


@app.after_request
def _no_cache(resp):
    resp.headers["Cache-Control"] = "no-store, max-age=0"
    resp.headers["Pragma"] = "no-cache"
    resp.headers["Expires"] = "0"
    return resp


@app.errorhandler(500)
def handle_500(e):
    return jsonify({"status": "error", "message": "Internal Server Error"}), 500


if __name__ == "__main__":
    print("Using data directory:", BASE_WORKING_DIR)

    if not os.path.exists(BASE_WORKING_DIR):
        print("WARNING: Base directory does not exist!")

    app.run(debug=True, host="0.0.0.0", port=PORT)
