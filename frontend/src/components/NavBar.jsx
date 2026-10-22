/**
 * @file NavBar.jsx
 * @summary Top navigation joining the original pyMSscreen workflow (prescreen -> extraction -> plotting)
 * with the HRMS NTS pages (non-target screening, semi-quant, product ions).
 */
import { NavLink } from 'react-router-dom';
import { useTheme } from './nts/ui';

const GROUPS = [
  ['Targeted prescreen (pyMSscreen)', [['/', 'Inputs'], ['/extraction', 'Extraction'], ['/plotting', 'Plotting']]],
  ['Non-target screening', [['/nts', 'NTS levels'], ['/semiquant', 'Semi-quant'], ['/product-ions', 'Product ions']]],
];

export default function NavBar() {
  const theme = useTheme();
  return (
    <nav style={{ display: 'flex', gap: 22, alignItems: 'center', flexWrap: 'wrap', padding: '8px 16px',
      background: theme.panel, borderBottom: `1px solid ${theme.border}`, position: 'sticky', top: 0, zIndex: 50 }}>
      <strong style={{ color: theme.text, fontSize: '1.05rem' }}>pyMSscreen-NTS</strong>
      {GROUPS.map(([group, links]) => (
        <div key={group} style={{ display: 'flex', gap: 6, alignItems: 'center' }}>
          <span style={{ color: theme.muted, fontSize: '0.75rem', textTransform: 'uppercase', letterSpacing: 0.5 }}>{group}</span>
          {links.map(([to, label]) => (
            <NavLink key={to} to={to} end
              style={({ isActive }) => ({ padding: '4px 10px', borderRadius: 6, fontWeight: 600, fontSize: '0.92rem', textDecoration: 'none',
                color: isActive ? '#fff' : theme.text, background: isActive ? theme.accent : 'transparent' })}>
              {label}
            </NavLink>
          ))}
        </div>
      ))}
    </nav>
  );
}
