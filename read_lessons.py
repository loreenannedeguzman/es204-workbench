from pathlib import Path
from pypdf import PdfReader
paths=[Path(r'F:\ES204\notes\ES204-Unit I.A-1_2026notes.pdf'),Path(r'F:\ES204\notes\ES204Unit II-A-2_2026notes.pdf'),Path(r'F:\ES204\notes\ES204Unit II-A-3_2026notes.pdf'),Path(r'C:\Users\Loreen Anne\Downloads\Tutorials Unit I.A-2 (Post-processing).pdf')]
out=Path(__file__).parent/'lesson_text';out.mkdir(exist_ok=True)
for i,path in enumerate(paths):
    reader=PdfReader(path)
    text='\n'.join(f'\n--- PAGE {j+1} ---\n'+(page.extract_text() or '') for j,page in enumerate(reader.pages))
    (out/f'{i+1}.txt').write_text(text,encoding='utf-8')
    print(path.name,len(reader.pages),'pages',len(text),'characters')
