import os
import base64
import json
import argparse
from datetime import datetime
from bs4 import BeautifulSoup

def build_inserts_dict(run_ids, barcodes, base_inserts_run, base_inserts_prefix):
    inserts = {}
    for run_id in run_ids:
        if run_id in base_inserts_run:
            inserts[run_id] = base_inserts_run[run_id]
    for run_id in run_ids:
        for bc in barcodes:
            prefix = f"{run_id}.{bc}"
            if prefix in base_inserts_prefix:
                inserts[prefix] = base_inserts_prefix[prefix]
    return inserts

def generate_html_report(prefixes, data_dir, output_path, sample_id, inserts):
    soup = BeautifulSoup("<html><head></head><body></body></html>", 'html.parser')

    style = soup.new_tag("style")
    style.string = """
    body { font-family: Arial, sans-serif; margin: 40px; }
    h1 { font-size: 2.2em; text-align: center; margin-bottom: 10px; }
    h2 { margin-top: 40px; border-bottom: 2px solid #ccc; padding-bottom: 5px; }
    table { border-collapse: collapse; margin-top: 15px; width: 100%; }
    th, td { border: 1px solid #ccc; padding: 8px; text-align: left; }
    tr:nth-child(even) { background-color: #f9f9f9; }
    img { display: block; margin: 20px auto; max-width: 600px; height: auto; border: 1px solid #ccc; padding: 4px; background-color: white; border-radius: 4px; }
    pre { background-color: #f4f4f4; padding: 10px; border-radius: 4px; overflow-x: auto; }
    .insert-block {
        background-color: #f0f8ff;
        border-left: 5px solid #339;
        padding: 12px 16px;
        margin: 20px 0;
        border-radius: 6px;
        font-size: 0.95em;
    }
    hr {
        margin: 30px auto;
        width: 60%;
        border: 1px solid #ccc;
    }
    .timestamp {
        text-align: center;
        font-size: 0.9em;
        color: #666;
        margin-bottom: 30px;
    }
    """
    soup.head.append(style)
    title_tag = soup.new_tag("title")
    title_tag.string = "Combined Report"
    soup.head.append(title_tag)
    h1 = soup.new_tag("h1")
    h1.string = f"Preliminary Processing Summary Report for Sample: {sample_id}"
    soup.body.append(h1)
    #soup.body.append(soup.new_tag("h1", string=f"📋 Preliminary Processing Summary Report for Sample: {sample_id}"))

    timestamp = soup.new_tag("div", **{"class": "timestamp"})
    timestamp.string = "Generated on: " + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    soup.body.append(timestamp)
    soup.body.append(soup.new_tag("hr"))

    def embed_image(path):
        with open(path, "rb") as f:
            return "data:image/png;base64," + base64.b64encode(f.read()).decode("utf-8")

    for prefix in prefixes:
        run_id, barcode_id = prefix.split('.', 1) if '.' in prefix else (prefix, None)
        section = soup.new_tag("div", style="margin-bottom: 3em;")

        # Section Title
        sec_title = soup.new_tag("h2")
        sec_title.string = f"Run: {run_id} | Barcode: {barcode_id}" if barcode_id else f"Run: {run_id}"
        section.append(sec_title)

        # Inserts
        for key, attr in [("before_txt", run_id), ("before_csv", prefix), ("before_png", prefix)]:
            html = inserts.get(attr, {}).get(key)
            if html:
                div = soup.new_tag("div", **{'class': 'insert-block'})
                div.append(BeautifulSoup(html, "html.parser"))
                section.append(div)

        # Text file
        txt_path = os.path.join(data_dir, f"{run_id}.ccs_report.txt")
        if os.path.exists(txt_path):
            with open(txt_path) as f:
                pre = soup.new_tag("pre")
                pre.string = f.read()
                section.append(pre)

        # CSV file
        csv_path = os.path.join(data_dir, f"{prefix}.skera.summary.csv")
        if os.path.exists(csv_path):
            with open(csv_path) as f:
                lines = f.readlines()
            table = soup.new_tag("table")
            for i, line in enumerate(lines):
                row = soup.new_tag("tr")
                for col in line.strip().split(","):
                    cell = soup.new_tag("th" if i == 0 else "td")
                    cell.string = col
                    row.append(cell)
                table.append(row)
            section.append(table)

        # PNGs
        for suffix in ["concat_hist.png", "ligations_heatmap.png", "readlen_hist.png"]:
            img_path = os.path.join(data_dir, f"{prefix}.{suffix}")
            if os.path.exists(img_path):
                img = soup.new_tag("img", src=embed_image(img_path))
                img['style'] = "max-width: 40%; display:block; margin-bottom: 1em; box-shadow: 2px 2px 6px rgba(0,0,0,0.15);"
                section.append(img)

        soup.body.append(section)

    # Extra Plot Section
    extra = soup.new_tag("div", style="margin-bottom: 3em;")
    extra.append(soup.new_tag("h2", string="Saturation Index - plotted using isoseq bcstats out"))
    extra.append(BeautifulSoup(f"<p>Saturation Index for <strong>{sample_id}</strong>.</p>", "html.parser"))
    extra.append(BeautifulSoup(f"<p>isoseq bcstats --method percentile --percentile 95 --json {sample_id}.bcstats.json -o {sample_id}.bcstats.tsv {sample_id}.merged.CB_sorted.bam </p>", "html.parser"))

    sample_png = os.path.join(data_dir, f"{sample_id}.saturation_index_plot.png")
    if os.path.exists(sample_png):
        with open(sample_png, "rb") as f:
            img_data = base64.b64encode(f.read()).decode("utf-8")
        img = soup.new_tag("img", src=f"data:image/png;base64,{img_data}")
        img['style'] = "max-width:60%;display:block;margin-bottom:1em;box-shadow: 2px 2px 6px rgba(0,0,0,0.15);"
        extra.append(img)
    soup.body.append(extra)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(str(soup))


def main():
    parser = argparse.ArgumentParser(description="Generate HTML QC report with embedded data and plots.")
    parser.add_argument("--prefixes", nargs='+', required=True, help="List of run_id.barcode_id prefixes (e.g., run1.bc01 run1.bc02)")
    parser.add_argument("--barcodes", nargs='+', required=True, help="List of barcode IDs (e.g., bc01 bc02)")
    parser.add_argument("--data-dir", required=True, help="Path to directory containing .txt, .csv, and .png files")
    parser.add_argument("--output", required=True, help="Output HTML path")
    parser.add_argument("--sample-id", required=True, help="Sample identifier")
    parser.add_argument("--base-inserts-run", default=None, help="Path to JSON file with base_inserts_run")
    parser.add_argument("--base-inserts-prefix", default=None, help="Path to JSON file with base_inserts_prefix")

    args = parser.parse_args()
    prefixes = args.prefixes
    barcodes = args.barcodes
    run_ids = sorted(set(p.split('.')[0] for p in prefixes))

    # Load insert dictionaries from JSON if provided
    base_inserts_run = json.load(open(args.base_inserts_run)) if args.base_inserts_run else {}
    base_inserts_prefix = json.load(open(args.base_inserts_prefix)) if args.base_inserts_prefix else {}

    inserts = build_inserts_dict(run_ids, barcodes, base_inserts_run, base_inserts_prefix)

    generate_html_report(
        prefixes=prefixes,
        data_dir=args.data_dir,
        output_path=args.output,
        sample_id=args.sample_id,
        inserts=inserts
    )


if __name__ == "__main__":
    main()

