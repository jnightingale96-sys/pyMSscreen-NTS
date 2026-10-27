import { BrowserRouter, Routes, Route } from 'react-router-dom';
import LandingPage from './components/LandingPage';
import ExtractionPage from './components/ExtractionPage';
import PlottingPage from './components/PlottingPage';
import NavBar from './components/NavBar';
import NtsPage from './components/NtsPage';
import SemiQuantPage from './components/SemiQuantPage';
import TargetedPage from './components/TargetedPage';
import ProductIonPage from './components/ProductIonPage';

export default function App() {
  return (
    <BrowserRouter>
      <NavBar />
      <Routes>
        <Route path="/" element={<LandingPage />} />
        <Route path="/extraction" element={<ExtractionPage />} />
        <Route path="/plotting" element={<PlottingPage />} />
        <Route path="/nts" element={<NtsPage />} />
        <Route path="/targeted" element={<TargetedPage />} />
        <Route path="/semiquant" element={<SemiQuantPage />} />
        <Route path="/product-ions" element={<ProductIonPage />} />
      </Routes>
    </BrowserRouter>
  );
}
