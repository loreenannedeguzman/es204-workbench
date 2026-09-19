from pathlib import Path
import pypdfium2 as pdfium
from PIL import Image,ImageDraw
paths=[r'F:\ES204\notes\ES204-Unit I.A-1_2026notes.pdf',r'F:\ES204\notes\ES204Unit II-A-2_2026notes.pdf',r'F:\ES204\notes\ES204Unit II-A-3_2026notes.pdf',r'C:\Users\Loreen Anne\Downloads\Tutorials Unit I.A-2 (Post-processing).pdf']
out=Path(__file__).parent/'lesson_images';out.mkdir(exist_ok=True)
for j,path in enumerate(paths,1):
    doc=pdfium.PdfDocument(path)
    for i in range(len(doc)):
        page=doc[i];image=page.render(scale=1.5).to_pil().convert('RGB')
        image.save(out/f'{j}-{i+1}.png')
    for start in range(0,len(doc),6):
        sheet=Image.new('RGB',(1500,2100),'#dddddd');draw=ImageDraw.Draw(sheet)
        for offset,i in enumerate(range(start,min(start+6,len(doc)))):
            im=Image.open(out/f'{j}-{i+1}.png');im.thumbnail((740,650))
            x=(offset%2)*750;y=(offset//2)*700
            sheet.paste(im,(x,y+25));draw.text((x+10,y+4),f'Document {j}, page {i+1}',fill='black')
        sheet.save(out/f'sheet-{j}-{start+1}.png')
    print(j,len(doc),flush=True)
