import yaml
import logging

from google.protobuf.json_format import ParseDict, SerializeToJsonError, ParseError
import ncbi.datasets.v1alpha1.reports.gene_pb2 as gene_report_pb2

logging.basicConfig(filename='process_markers.log', level=logging.DEBUG)
logger = logging.getLogger(__name__)


class DatasetsReportReader():
    # Load yaml in 'buf' and parse resulting dictionary into protobuf 'schema_pb'
    def _load_and_parse_report(self, buf, schema_pb, from_array=False):
        try:
            report_dict = yaml.safe_load(buf)
            if not report_dict:
                logger.error("Empty report from file")
            try:
                if from_array:
                    ParseDict(report_dict[0], schema_pb, ignore_unknown_fields=False)
                else:
                    ParseDict(report_dict, schema_pb, ignore_unknown_fields=False)
            except (SerializeToJsonError, ParseError) as e:
                logger.error(f"Error converting yaml to schema: {e}")
        except yaml.YAMLError as e:
            logger.error(f"Error while loading yaml: {e}")
        return schema_pb

    # return full gene report
    def gene_report_from_file(self, report_file_name):
        schema_pb = gene_report_pb2.GeneDescriptors()
        with open(report_file_name) as fh:
            self._load_and_parse_report(fh.read(), schema_pb)
        return schema_pb
