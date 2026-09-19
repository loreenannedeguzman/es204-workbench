import json
import unittest
from numerics import execute
from final_values import meets, precision_warning


class FinalValueChecks(unittest.TestCase):
    def inputs(self,**changes):
        p=dict(formula='20*cos(x)-40*y',x0='0',xf='1',y0='0',
               grid_mode='Powers of ten n=10^p',first_p='1',last_p='6',
               output_mode='Final values only',methods=['Euler','Predictor-Corrector','RK3'],
               iterations='1',exact='20/1601*(40*cos(x)+sin(x)-40*exp(-40*x))',reference='',
               stop_metric='Final true error %',stop_tolerance='1e-5',comparison='<')
        p.update(changes)
        return p

    def run_case(self,**changes):
        return execute({'action':'solve','inputs':self.inputs(**changes)})

    def test_million_step_and_separate_stops(self):
        result=self.run_case()
        rows=result['tables'][0]['rows']
        passed=[r for r in rows if r[7]]
        self.assertEqual([r[1] for r in passed],[1000000,10000,1000])
        expected=[.280494234982769,.280494230629948,.280494231843986]
        for row,value in zip(passed,expected):self.assertAlmostEqual(row[3],value,places=13)
        self.assertLess(passed[0][4],0);self.assertGreater(passed[1][4],0)
        self.assertEqual(len(rows),13)
        self.assertEqual(result['datasets'],[])
        self.assertLess(len(json.dumps(result,allow_nan=False)),12000)

    def test_rounded_reference_retained(self):
        result=self.run_case(exact='',reference='1, 0.280494',last_p='4')
        self.assertTrue(result['warnings'])
        self.assertFalse(any(row[7] for row in result['tables'][0]['rows']))
        for row in result['tables'][0]['rows']:
            self.assertAlmostEqual(row[4],(.280494-row[3])/.280494*100)

    def test_explicit_table_over_formula(self):
        result=self.run_case(reference='1,0.280494',final_reference_source='Given final reference table',methods=['RK3'],last_p='3')
        self.assertTrue(result['warnings']);self.assertFalse(result['tables'][0]['rows'][-1][7])

    def test_exact_equality_strict_and_inclusive(self):
        changes=dict(formula='0',x0='0',xf='1',y0='1',exact='',reference='1,2',
                     methods=['Euler'],stop_tolerance='50',first_p='0',last_p='1')
        strict=self.run_case(**changes,comparison='<')
        inclusive=self.run_case(**changes,comparison='<=')
        self.assertEqual(len(strict['tables'][0]['rows']),2)
        self.assertFalse(strict['tables'][0]['rows'][-1][7])
        self.assertEqual(len(inclusive['tables'][0]['rows']),1)
        self.assertTrue(inclusive['tables'][0]['rows'][0][7])

    def test_final_modes_and_original_results(self):
        for mode,steps in [('Step sizes h','.1,.05'),('Interval counts n','10,20'),('Refinement h0 / 2^m','.1')]:
            changes=dict(grid_mode=mode,steps=steps,max_m='1',last_p='2',methods=['RK3'],stop_metric='No early stop')
            compact=self.run_case(**changes)
            full=self.run_case(**changes,output_mode='Full tables')
            for row,data in zip(compact['tables'][0]['rows'],full['datasets']):
                self.assertAlmostEqual(row[3],data['y'][-1],places=13)
            self.assertTrue(any('signed true error %' in t['columns'] for t in full['tables']))

    def test_missing_zero_reference_and_limits(self):
        for changes in [dict(exact='',reference=''),dict(exact='',reference='1,0'),dict(last_p='8')]:
            with self.assertRaises(ValueError):self.run_case(**changes)
        self.assertIsNone(precision_warning('0.280494231349575',1e-5))
        self.assertFalse(meets(None,1,'<'))


if __name__=='__main__':unittest.main()
