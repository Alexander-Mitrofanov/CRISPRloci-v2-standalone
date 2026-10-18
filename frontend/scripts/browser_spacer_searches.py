#!/usr/bin/env python3
"""Browser acceptance for both spacer-search directions using synthetic evidence."""
import argparse
import copy
import json
from pathlib import Path
from urllib.parse import urlsplit
from playwright.sync_api import expect, sync_playwright
from browser_protospacer import check as check_forward

ROOT = Path(__file__).resolve().parents[1]
TOKEN = 'a' * 43
FASTA = '>spacer_A\nACGTTGCAACGATTCGAGTC\n'
TABLE = b'query_id\tsubject_accession\nspacer_A\tsynthetic_A\n'


def check(origin, output):
    check_forward(origin, output / 'forward')
    output.mkdir(parents=True, exist_ok=True)
    fixture = json.loads((ROOT / 'tests/fixtures/viral.json').read_text())
    fixture['artifacts'] = [{'artifact_id':'viral-table','filename':'viral-matches.tsv', 'label':'Complete viral matches','media_type':'text/tab-separated-values','size_bytes':len(TABLE)}]
    state = {'job':fixture,'submissions':[], 'modes':['protospacer','viral_search']}
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        context = browser.new_context(viewport={'width':1440,'height':1100}, accept_downloads=True)
        page = context.new_page()
        page.on('pageerror', lambda e: errors.append(str(e)))
        page.add_init_script("window.cspViolations=[];document.addEventListener('securitypolicyviolation',e=>window.cspViolations.push(e.violatedDirective))")
        def api(route):
            req = route.request
            path = urlsplit(req.url).path
            headers = {'Access-Control-Allow-Origin':origin,'Access-Control-Allow-Headers':'authorization,content-type'}
            if req.method == 'OPTIONS':
                route.fulfill(status=204,headers=headers)
                return
            if path.endswith('/health'):
                data = {'version':'2.3.0'}
            elif path.endswith('/config'):
                data = {'modes':state['modes'],'max_total_bases':10000,'max_record_bases':10000,'max_records':20,'max_request_bytes':20000}
            elif path.endswith('/jobs') and req.method == 'POST':
                body = req.post_data_json
                assert body['mode']=='viral_search' and body['viral_max_mismatches']==0 and body['sequence']==FASTA
                assert not {'reference','reference_path','blastn'} & body.keys()
                state['submissions'].append(body)
                data = {'job_id':fixture['job_id'],'access_token':TOKEN,'job':state['job']}
            elif '/jobs/' in path:
                assert req.headers.get('authorization')=='Bearer '+TOKEN and TOKEN not in req.url
                if path.endswith('/artifacts/viral-table'):
                    route.fulfill(status=200,headers=headers,content_type='text/tab-separated-values',body=TABLE)
                    return
                data = state['job']
            else:
                raise AssertionError(path)
            route.fulfill(status=202 if req.method=='POST' else 200,headers=headers,content_type='application/json',body=json.dumps(data))
        page.route('https://analysis.example.org/**',api)
        page.goto(origin)
        page.get_by_role('link',name='Spacer searches',exact=True).click()
        page.get_by_label('Contigs or small genomes',exact=True).fill('>my_contig\nACGTACGT\n')
        page.get_by_role('button',name='Find viruses for my spacers',exact=True).click()
        expect(page.get_by_role('heading',name='Spacer searches',exact=True)).to_be_visible()
        page.get_by_label('Spacer sequences',exact=True).fill(FASTA)
        page.get_by_role('button',name='Find spacers in my sequence',exact=True).click()
        expect(page.get_by_label('Contigs or small genomes',exact=True)).to_have_value('>my_contig\nACGTACGT\n')
        page.get_by_role('button',name='Find viruses for my spacers',exact=True).click()
        expect(page.get_by_label('Spacer sequences',exact=True)).to_have_value(FASTA)
        page.get_by_label('Maximum substitutions',exact=False).select_option('0')
        page.evaluate('window.scrollTo(0,0)')
        page.screenshot(path=output/'search-desktop.png',full_page=True)
        page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth<=window.innerWidth+1')
        page.evaluate('window.scrollTo(0,0)')
        page.screenshot(path=output/'search-mobile.png',full_page=True)
        page.set_viewport_size({'width':1440,'height':1100})
        page.get_by_role('button',name='Search references',exact=True).click()
        overview = page.get_by_role('tabpanel', name='Overview', exact=True)
        expect(overview.get_by_text('BLAST retrieval is heuristic.', exact=False)).to_be_visible()
        expect(overview.get_by_role('table')).to_have_count(0)
        expect(overview.get_by_label('Filter displayed matches by spacer, accession or title', exact=True)).to_have_count(0)
        page.get_by_role('tab',name='Viral matches',exact=True).click()
        panel=page.get_by_role('tabpanel',name='Viral matches',exact=True)
        expect(panel.get_by_role('row')).to_have_count(4)
        page.evaluate('window.scrollTo(0,0)')
        page.screenshot(path=output/'viral-desktop.png',full_page=True)
        panel.get_by_label('Filter displayed matches by spacer, accession or title',exact=True).fill('no-such-virus')
        expect(panel.get_by_text('No displayed matches contain that text.',exact=True)).to_be_visible()
        page.reload()
        page.get_by_role('tab',name='Viral matches',exact=True).click()
        expect(panel.get_by_role('row')).to_have_count(4)
        page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth<=window.innerWidth+1')
        page.evaluate('window.scrollTo(0,0)')
        page.screenshot(path=output/'viral-mobile.png',full_page=True)
        page.get_by_role('tab',name='Files & methods',exact=True).click()
        with page.expect_download() as result:
            page.get_by_role('button',name='Complete viral matches',exact=False).click()
        result.value.save_as(output/'viral-matches.tsv')
        assert (output/'viral-matches.tsv').read_bytes()==TABLE
        for evidence, text in [(None,'A completed viral search was not reported.'),({'available':True,'execution_completed':True,'outcome':'no_eligible_queries'},'No eligible spacers to search.'),({'available':True,'execution_completed':True,'outcome':'no_raw_hsps'},'No full-length matches passed this search policy.')]:
            state['job']=copy.deepcopy(fixture)
            state['job']['summary']['viral_search']=evidence
            page.reload()
            expect(page.get_by_role('tabpanel',name='Overview',exact=True).get_by_text(text,exact=False)).to_be_visible()
        assert len(state['submissions'])==1 and not errors and not page.evaluate('window.cspViolations')
        (output/'acceptance.json').write_text(json.dumps({'status':'passed','api':'synthetic intercepted fixtures','checks':['both_directions','independent_drafts','zero_substitution_payload','mobile_forms_results','candidate_evidence','filter','recovery_reload','authenticated_download','missing_skipped_negative_states','no_console_or_csp_errors']},indent=2)+'\n')
        browser.close()
        print('Bidirectional spacer browser acceptance passed:',output)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--origin',default='http://127.0.0.1:4187')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    check(args.origin.rstrip('/'),args.output)
