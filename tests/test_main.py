import logging
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from main import configure_logging, filter_fastq, parse_args


FASTQ_TEMPLATE = '@{name}\n{seq}\n+\n{qual}\n'


def write_fastq(path, records):
    path.write_text(''.join(
        FASTQ_TEMPLATE.format(name=name, seq=seq, qual=qual)
        for name, seq, qual in records
    ), encoding='utf-8')


def make_logger(path):
    logging.getLogger('biogen').handlers.clear()
    return configure_logging(str(path))


class TestFilterFastq:
    def test_filter_fastq_creates_output_file(self, tmp_path):
        input_path = tmp_path / 'input.fastq'
        write_fastq(input_path, [('read1', 'ACGT', 'IIII')])

        output_path = filter_fastq(
            str(input_path),
            'out.fastq',
            logger=make_logger(tmp_path / 'biogen.log'),
        )

        assert output_path is not None
        assert Path(output_path).exists()
        assert 'filtered' in output_path

    def test_filter_fastq_gc_bounds_filters_records(self, tmp_path):
        input_path = tmp_path / 'input.fastq'
        write_fastq(input_path, [
            ('high_gc', 'GCGC', 'IIII'),
            ('low_gc', 'ATAT', 'IIII'),
        ])

        output_path = filter_fastq(
            str(input_path),
            'out.fastq',
            gc_bounds=(75, 100),
            logger=make_logger(tmp_path / 'biogen.log'),
        )

        text = Path(output_path).read_text(encoding='utf-8')
        assert '@high_gc' in text
        assert '@low_gc' not in text

    def test_filter_fastq_length_bounds_filters_records(self, tmp_path):
        input_path = tmp_path / 'input.fastq'
        write_fastq(input_path, [
            ('short', 'AT', 'II'),
            ('long', 'ATATAT', 'IIIIII'),
        ])

        output_path = filter_fastq(
            str(input_path),
            'out.fastq',
            length_bounds=(0, 3),
            logger=make_logger(tmp_path / 'biogen.log'),
        )

        text = Path(output_path).read_text(encoding='utf-8')
        assert '@short' in text
        assert '@long' not in text

    def test_filter_fastq_quality_threshold_filters_records(self, tmp_path):
        input_path = tmp_path / 'input.fastq'
        write_fastq(input_path, [
            ('good', 'ATGC', 'IIII'),
            ('bad', 'ATGC', '!!!!'),
        ])

        output_path = filter_fastq(
            str(input_path),
            'out.fastq',
            quality_threshold=20,
            logger=make_logger(tmp_path / 'biogen.log'),
        )

        text = Path(output_path).read_text(encoding='utf-8')
        assert '@good' in text
        assert '@bad' not in text

    def test_filter_fastq_overwrite_disabled_returns_none(self, tmp_path):
        input_path = tmp_path / 'input.fastq'
        write_fastq(input_path, [('read1', 'ACGT', 'IIII')])

        filtered_dir = tmp_path / 'filtered'
        filtered_dir.mkdir()
        existing_output = filtered_dir / 'out.fastq'
        existing_output.write_text('existing', encoding='utf-8')

        result = filter_fastq(
            str(input_path),
            'out.fastq',
            overwrite=False,
            logger=make_logger(tmp_path / 'biogen.log'),
        )

        assert result is None


class TestCliAndLogging:
    def test_parse_args_returns_expected_values(self):
        args = parse_args([
            'input.fastq',
            'out.fastq',
            '--gc-min', '10',
            '--gc-max', '90',
            '--length-min', '5',
            '--length-max', '50',
            '--quality-threshold', '25',
            '--overwrite',
            '--log-file', 'app.log',
        ])

        assert args.input_fastq == 'input.fastq'
        assert args.output_fastq == 'out.fastq'
        assert args.gc_min == 10.0
        assert args.gc_max == 90.0
        assert args.length_min == 5
        assert args.length_max == 50
        assert args.quality_threshold == 25
        assert args.overwrite is True
        assert args.log_file == 'app.log'

    def test_filter_fastq_logs_error_for_missing_input(self, tmp_path):
        log_path = tmp_path / 'biogen.log'
        logger = make_logger(log_path)

        result = filter_fastq(
            str(tmp_path / 'missing.fastq'),
            'out.fastq',
            logger=logger,
        )

        assert result is None
        content = log_path.read_text(encoding='utf-8')
        assert 'Input FASTQ file not found' in content

    def test_filter_fastq_writes_info_to_log_file(self, tmp_path):
        log_path = tmp_path / 'biogen.log'
        logger = make_logger(log_path)

        input_path = tmp_path / 'input.fastq'
        write_fastq(input_path, [('read1', 'ACGT', 'IIII')])

        output_path = filter_fastq(
            str(input_path),
            'out.fastq',
            logger=logger,
        )

        assert output_path is not None
        content = log_path.read_text(encoding='utf-8')
        assert 'Starting FASTQ filter' in content
        assert 'Filtered FASTQ written to' in content
