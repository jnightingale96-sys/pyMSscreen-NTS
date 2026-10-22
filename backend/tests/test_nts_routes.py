"""Routes that connect the web app to the HRMS NTS engine (skipped when engine/ is not installed, e.g. Docker CI)."""
import pytest

pytest.importorskip("hrms_nts")

from web_app import app  # noqa: E402


@pytest.fixture
def client():
    app.config["TESTING"] = True
    return app.test_client()


@pytest.fixture
def run_dir(tmp_path):
    d = tmp_path / "NTS_out_test"
    d.mkdir()
    (d / "level_summary.csv").write_text("sample,1,2a,2b,3,4,5\nS1,0,1,0,0,1,0\n")
    (d / "compound_summary.csv").write_text(
        "name,ik1,confidence_score,area,max_int,rt_apex,ppm_error,n_samples_detected,formula,smiles,best_level\n"
        "Carbamazepine,FFGPTBGBLSHEPO,90,2.0e7,5e6,8.0,0.5,1,C15H12N2O,NC(=O)N1c2ccccc2C=Cc2ccccc12,2a\n"
        "Caffeine,RYYVLZVUVIJVGH,40,1.0e6,2e5,3.0,1.2,1,C8H10N4O2,Cn1cnc2c1c(=O)n(C)c(=O)n2C,4\n")
    (d / "top_candidates_per_feature.csv").write_text(
        "sample,feature_id,polarity,level,name,smiles,rt_apex,area,max_int\n"
        "S1,S1_F000001,positive,2a,Carbamazepine,NC(=O)N1c2ccccc2C=Cc2ccccc12,8.0,2.0e7,5e6\n"
        "S1,S1_F000002,positive,4,Caffeine,Cn1cnc2c1c(=O)n(C)c(=O)n2C,3.0,1.0e6,2e5\n")
    return d


def test_info(client):
    r = client.get("/nts/info").get_json()
    assert r["version"] and "water" in r["matrices"] and r["checks"]


def test_runs_and_results(client, run_dir):
    runs = client.get("/nts/runs", query_string={"folder": str(run_dir.parent)}).get_json()["runs"]
    assert [x["name"] for x in runs] == ["NTS_out_test"]
    r = client.get("/nts/results", query_string={"run_dir": str(run_dir), "table": "compounds", "levels": "1,2a"}).get_json()
    assert r["total"] == 1 and r["rows"][0]["name"] == "Carbamazepine"


def test_file_route_stays_inside_results_folder(client, run_dir):
    assert client.get("/nts/file", query_string={"run_dir": str(run_dir), "name": "level_summary.csv"}).status_code == 200
    (run_dir.parent / "secret.csv").write_text("x\n")
    assert client.get("/nts/file", query_string={"run_dir": str(run_dir), "name": "../secret.csv"}).status_code == 404


def test_semiquant_single_reference(client, run_dir):
    r = client.post("/nts/semiquant", json={"run_dir": str(run_dir), "mode": "single", "matrix": "extract", "levels": ["2a", "4"],
                                            "single_reference": {"name": "Ref", "conc_ug_mL": 1.0, "area": 1.0e7}}).get_json()
    conc = {row["name"]: row["conc_extract"] for row in r["rows"]}
    assert conc["Carbamazepine"] == pytest.approx(2.0)          # 1 ug/mL x 2e7 / 1e7
    assert conc["Caffeine"] == pytest.approx(0.1)
    assert (run_dir / "semiquant_single_extract.csv").exists()


def test_export_to_prescreen(client, run_dir):
    r = client.post("/nts/export_to_prescreen", json={"run_dir": str(run_dir), "levels": ["2a"]}).get_json()
    assert r["n"] == 1
    text = open(r["path"], encoding="utf-8").read()
    assert text.startswith("ID,Name,SMILES,mz") and "Carbamazepine" in text


def test_predict_ions_and_structure(client):
    r = client.post("/nts/predict_ions", json={"smiles": "NC(=O)N1c2ccccc2C=Cc2ccccc12", "adducts": ["[M+H]+"],
                                               "experimental_peaks": "194.0964 100\n192.0808 30"}).get_json()
    res = r["results"][0]
    assert r["formula"] == "C15H12N2O" and abs(res["precursor_mz"] - 237.1022) < 1e-3
    assert res["ions"] and res["transitions"] and "bi_score" in res["experimental"]
    assert client.post("/nts/predict_ions", json={"smiles": "not a smiles"}).status_code == 400
    svg = client.get("/nts/structure.svg", query_string={"smiles": "CCO"})
    assert svg.status_code == 200 and b"<svg" in svg.data
