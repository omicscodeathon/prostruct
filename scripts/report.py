import os
import yagmail
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, TableStyle, Spacer
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.charts.legends import Legend
from reportlab.graphics.widgets.markers import makeMarker
from reportlab.graphics import renderPDF

from analyser import basic_analysis

#
# def generate_report(filename,job_name,sequence_analysis_result):
#     document = SimpleDocTemplate(filename=filename,pagesize=A4)
#     elements = []
#
#     styles = getSampleStyleSheet()
#     title = Paragraph(f"{job_name} Report",style=styles['Title'])
#     introduction = Paragraph(f"Analysis for {job_name}",style=styles["BodyText"])
#     elements.append(title)
#     elements.append(Spacer(1, 12))
#     elements.append(introduction)
#     elements.append(Spacer(1, 12))
#
#     keys = list(sequence_analysis_result.keys())
#     count_talbe = [["Amino Acid",'Count']]
#     for i in range(len(sequence_analysis_result)):
#         if keys[i] == "amino acid count":
#             count_data = sequence_analysis_result[keys[i]]
#             for key,value in count_data.items():
#                 count_talbe.append([key,value])
#             table = Table(count_talbe,colWidths=[210,210])
#             table.setStyle(TableStyle([
#                     ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
#                     ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
#                     ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
#                     ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
#                     ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
#                     ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
#                     ('GRID', (0, 0), (-1, -1), 1, colors.black)
#                 ]))
#             elements.append(table)
#
#     for i in range(len(sequence_analysis_result)):
#         if keys[i] != "amino acid count":
#             data = f"{keys[i]}: {sequence_analysis_result[keys[i]]}"
#             sub_element = Paragraph(data,style=styles["BodyText"])
#             elements.append(Spacer(1,2))
#             elements.append(sub_element)
#
#
#     document.build(elements)
#
    # return report



def send_report(recipient,job_name,modeled_structures_directory="modeled_structures",report=None):
    attachments = []
    receiver = recipient
    body = """
            Dear Prostruc user,
            
            Kindly find attached your predicted structures and a comprehensive report.
            Please feel free to reach out to us if the need arises.
            
            Best,
            Prostruc Team
    """
    try:
        if os.path.exists(modeled_structures_directory):
            for path,directories,files in os.walk(modeled_structures_directory):
                for target_model in files:
                    target_model_filepath = os.path.join(modeled_structures_directory,target_model)
                    attachments.append(target_model_filepath)



        yag = yagmail.SMTP("prostruct.omic@gmail.com",password="ztzwfcurbjfkjxdj")
        yag.send(to=receiver, subject=f"Prostruc Job Result :: {job_name}", contents=body, attachments=attachments)
        print("[*] Report sent")

    except Exception as error:
        print(error)
        return False


# analysis_result = basic_analysis(sequence="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN")
# generate_report(filename="pro2.pdf",job_name="pro2",sequence_analysis_result=analysis_result)
# print(analysis_result.keys())
