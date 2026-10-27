import axios from 'axios';

// Docker: separate backend container on :8080. Native Windows build (`npm run build:windows`): Flask serves
// this app and the API from the same origin.
export const API_BASE =
  import.meta.env.VITE_API_BASE ||
  (import.meta.env.MODE === 'windows' ? window.location.origin : 'http://localhost:8080');

// ---- HRMS NTS engine (non-target screening, semi-quant, product ions) ----
export const ntsInfo = () => axios.get(`${API_BASE}/nts/info`);
export const ntsListMzml = (folder) => axios.get(`${API_BASE}/nts/list_mzml`, { params: { folder } });
export const ntsStartRun = (payload) => axios.post(`${API_BASE}/nts/start_run`, payload);
export const ntsStartCalibrate = (payload) => axios.post(`${API_BASE}/nts/start_calibrate`, payload);
export const ntsStartTargeted = (payload) => axios.post(`${API_BASE}/nts/start_targeted`, payload);
export const ntsTargetedResults = (run_dir) => axios.get(`${API_BASE}/nts/targeted_results`, { params: { run_dir } });
export const ntsJob = (job_id) => axios.get(`${API_BASE}/nts/job`, { params: { job_id } });
export const ntsCancel = (job_id) => axios.post(`${API_BASE}/nts/cancel`, { job_id });
export const ntsRuns = (folder) => axios.get(`${API_BASE}/nts/runs`, { params: { folder } });
export const ntsResults = (params) => axios.get(`${API_BASE}/nts/results`, { params });
export const ntsPlots = (run_dir) => axios.get(`${API_BASE}/nts/plots`, { params: { run_dir } });
export const ntsFileUrl = (run_dir, name, download = false) =>
  `${API_BASE}/nts/file?${new URLSearchParams({ run_dir, name, download: download ? '1' : '0' })}`;
export const ntsExportToPrescreen = (payload) => axios.post(`${API_BASE}/nts/export_to_prescreen`, payload);
export const ntsGetCalibrants = () => axios.get(`${API_BASE}/nts/calibrants`);
export const ntsSaveCalibrants = (rows) => axios.post(`${API_BASE}/nts/calibrants`, { rows });
export const ntsSemiquant = (payload) => axios.post(`${API_BASE}/nts/semiquant`, payload);
export const ntsPredictIons = (payload) => axios.post(`${API_BASE}/nts/predict_ions`, payload);
export const ntsStructureUrl = (smiles, size = 260) =>
  `${API_BASE}/nts/structure.svg?${new URLSearchParams({ smiles, size })}`;
export const apiError = (err) => err?.response?.data?.error || err?.response?.data?.message || err?.message || 'Unknown error';

export const saveState = (state) => {
  return axios.post(`${API_BASE}/save_state`, state);
};

export const loadState = (payload) => {
  return axios.post(`${API_BASE}/load_state`, payload);
};

export const generateTable = (payload) => {
  return axios.post(`${API_BASE}/generate_table`, payload, {
    headers: {
      'Content-Type': 'application/json',
    }
  });
};

export const saveExtractConfig = async (config) => {
  return axios.post(`${API_BASE}/save_extract_config`, config);
};

export const loadExtractionConfig = async (workingDir) => {
  try {
    const params = new URLSearchParams({ working_directory: workingDir });
    const res = await axios.get(`${API_BASE}/load_extraction_config`, { params });
    return parseFloat(res.data.ret_time_shift_tol) || 0.5;
  } catch {
    return 0.5;
  }
};

export const exportSummaryPDF = async (compounds, working_directory) => {
  // Properly structure payload for backend
  const payload = {
    working_directory,
    compound_groups: compounds.reduce((groups, compound) => {
      const group = groups.find(g => g.group_name === compound.ID);
      if (group) {
        group.compounds.push(compound);
      } else {
        groups.push({
          group_name: compound.ID,
          compounds: [compound]
        });
      }
      return groups;
    }, [])
  };

  try {
    const response = await axios.post(`${API_BASE}/export_summary_pdf`, payload, {
      responseType: 'blob'
    });
    return response.data;
  } catch (error) {
    console.error('PDF export failed:', error);
    throw new Error('PDF export failed');
  }
};

export const saveQAFlags = (working_directory, flags) => {
  return axios.post(`${API_BASE}/save_qa_flags`, {
    working_directory,
    flags
  });
};
