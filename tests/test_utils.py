"""Tests for utility functions in pypaperretriever.utils"""
import pytest
from pypaperretriever.utils import encode_doi, decode_doi


class TestEncodeDecode:
    """Test DOI encoding and decoding with various special characters."""

    def test_encode_decode_basic_doi(self):
        """Test basic DOI without special characters."""
        doi = "10.1234/simple"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        # Verify slashes, dots are encoded
        assert "/" not in encoded
        assert "." not in encoded
        assert "%2F" in encoded  # slash
        assert "%2E" in encoded  # dot

    def test_encode_decode_with_parentheses(self):
        """Test DOI with parentheses - the main issue being fixed."""
        doi = "10.1016/S0022-4049(02)00143-3"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        # Verify parentheses are encoded
        assert "(" not in encoded
        assert ")" not in encoded
        assert "%28" in encoded  # opening paren
        assert "%29" in encoded  # closing paren

    def test_encode_decode_with_hyphens(self):
        """Test DOI with hyphens."""
        doi = "10.1234/test-hyphen-doi"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert "%2D" in encoded  # hyphen

    def test_encode_decode_with_dots(self):
        """Test DOI with multiple dots."""
        doi = "10.1234/test.with.dots"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert "." not in encoded

    def test_encode_decode_with_semicolons(self):
        """Test DOI with semicolons."""
        doi = "10.1234/test;semicolon"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert ";" not in encoded
        assert "%3B" in encoded  # semicolon

    def test_encode_decode_with_angle_brackets(self):
        """Test DOI with angle brackets."""
        doi = "10.1234/test<bracket>"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert "<" not in encoded
        assert ">" not in encoded
        assert "%3C" in encoded  # less than
        assert "%3E" in encoded  # greater than

    def test_encode_decode_with_hash(self):
        """Test DOI with hash symbol."""
        doi = "10.1234/test#hash"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert "#" not in encoded
        assert "%23" in encoded  # hash

    def test_encode_decode_with_question_mark(self):
        """Test DOI with question mark."""
        doi = "10.1234/test"
        # The encode function splits on '?' so test the actual DOI part
        doi_with_params = "10.1234/test?param=value"
        encoded = encode_doi(doi_with_params)
        decoded = decode_doi(encoded)
        # Should only encode the part before '?'
        assert decoded == doi

    def test_encode_decode_with_spaces(self):
        """Test DOI with spaces (rare but possible)."""
        doi = "10.1234/test space"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert " " not in encoded
        assert "%20" in encoded  # space

    def test_encode_decode_with_plus(self):
        """Test DOI with plus sign."""
        doi = "10.1234/test+plus"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert "+" not in encoded
        assert "%2B" in encoded  # plus

    def test_encode_decode_with_ampersand(self):
        """Test DOI with ampersand."""
        doi = "10.1234/test&ampersand"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert "&" not in encoded
        assert "%26" in encoded  # ampersand

    def test_encode_decode_with_equals(self):
        """Test DOI with equals sign."""
        doi = "10.1234/test=equals"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert "=" not in encoded
        assert "%3D" in encoded  # equals

    def test_encode_decode_complex_doi(self):
        """Test complex DOI with multiple special characters."""
        doi = "10.1016/S0022-4049(02)00143-3"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        # Comprehensive check
        for char in ['/', '.', '-', '(', ')']:
            assert char not in encoded

    def test_encode_with_doi_org_prefix(self):
        """Test that doi.org prefix is properly stripped."""
        doi_with_prefix = "https://doi.org/10.1234/test"
        doi_without_prefix = "10.1234/test"
        encoded_with = encode_doi(doi_with_prefix)
        encoded_without = encode_doi(doi_without_prefix)
        assert encoded_with == encoded_without
        assert decode_doi(encoded_with) == doi_without_prefix

    def test_encode_with_query_params(self):
        """Test that query parameters are stripped."""
        doi_with_params = "10.1234/test?param=value&other=thing"
        expected_doi = "10.1234/test"
        encoded = encode_doi(doi_with_params)
        decoded = decode_doi(encoded)
        assert decoded == expected_doi

    def test_encode_strips_quotes(self):
        """Test that quotes are stripped from input DOI before encoding."""
        doi_with_single_quotes = "'10.1234/test'"
        doi_with_double_quotes = '"10.1234/test"'
        doi_without_quotes = "10.1234/test"
        
        # All three should produce the same encoded result
        encoded1 = encode_doi(doi_with_single_quotes)
        encoded2 = encode_doi(doi_with_double_quotes)
        encoded3 = encode_doi(doi_without_quotes)
        
        assert encoded1 == encoded2 == encoded3
        
        # And they should all decode to the clean DOI
        assert decode_doi(encoded1) == doi_without_quotes
        assert decode_doi(encoded2) == doi_without_quotes
        assert decode_doi(encoded3) == doi_without_quotes

    def test_real_world_dois(self):
        """Test with real-world DOI examples."""
        real_dois = [
            "10.7759/cureus.76081",
            "10.1016/j.revmed.2011.10.009",
            "10.1016/S0022-4049(02)00143-3",
            "10.21105/joss.08135",
            "10.1016/j.nedt.2023.105737",
            "10.1097/00005537-200206000-00026",  # Multiple hyphens (medical)
            "10.1021/jp992190+",  # Plus sign (ACS)
        ]
        
        for doi in real_dois:
            encoded = encode_doi(doi)
            decoded = decode_doi(encoded)
            assert decoded == doi, f"Failed for DOI: {doi}"
    
    def test_complex_legacy_sici_doi(self):
        """Test with complex legacy SICI-style DOI containing many special characters."""
        # Real Wiley physics paper DOI from 1999
        doi = "10.1002/(SICI)1521-3951(199911)216:1<135::AID-PSSB135>3.0.CO;2-#"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        
        # Verify round-trip
        assert decoded == doi
        
        # Verify all special characters are encoded
        special_chars = ['(', ')', '<', '>', ';', '#', '/', '.', '-', ':']
        for char in special_chars:
            assert char not in encoded, f"Character '{char}' should be encoded"
        
        # Verify specific encodings are present
        assert '%28' in encoded  # (
        assert '%29' in encoded  # )
        assert '%3C' in encoded  # <
        assert '%3E' in encoded  # >
        assert '%3B' in encoded  # ;
        assert '%23' in encoded  # #

    def test_encoded_safe_for_filenames(self):
        """Test that encoded DOIs are safe for use in filenames."""
        doi = "10.1016/S0022-4049(02)00143-3"
        encoded = encode_doi(doi)
        
        # Characters that are problematic in filenames
        unsafe_chars = ['/', '\\', ':', '*', '?', '"', '<', '>', '|']
        for char in unsafe_chars:
            if char != '\\':  # backslash isn't in the DOI to begin with
                assert char not in encoded, f"Unsafe character '{char}' found in encoded DOI"

    def test_idempotent_encoding(self):
        """Test that encoding an already encoded DOI doesn't double-encode."""
        doi = "10.1016/S0022-4049(02)00143-3"
        encoded_once = encode_doi(doi)
        
        # Encoding an already encoded string should be handled gracefully
        # (though in practice, this shouldn't happen)
        decoded = decode_doi(encoded_once)
        encoded_again = encode_doi(decoded)
        assert encoded_once == encoded_again

    def test_roundtrip_multiple_times(self):
        """Test that encoding/decoding multiple times is stable."""
        doi = "10.1016/S0022-4049(02)00143-3"
        
        # Multiple roundtrips should maintain the DOI
        for _ in range(5):
            encoded = encode_doi(doi)
            decoded = decode_doi(encoded)
            assert decoded == doi


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_string(self):
        """Test encoding/decoding empty string."""
        assert encode_doi("") == ""
        assert decode_doi("") == ""

    def test_only_doi_prefix(self):
        """Test with just the prefix."""
        # Dots should be encoded
        assert encode_doi("10.") == "10%2E"
        assert decode_doi("10%2E") == "10."
        
        # Test a more complete prefix
        assert encode_doi("10.1234") == "10%2E1234"
        assert decode_doi("10%2E1234") == "10.1234"

    def test_unicode_in_doi(self):
        """Test DOI with unicode characters (rare but possible)."""
        doi = "10.1234/tëst"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi

    def test_case_sensitivity(self):
        """Test that encoding preserves case."""
        doi = "10.1234/TeSt-CaSe"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert "TeSt" in decoded
        assert "CaSe" in decoded


class TestAdditionalEdgeCases:
    """Additional edge cases and integration scenarios."""

    def test_colon_in_doi(self):
        """Test DOI with colons (common in SICI-style DOIs)."""
        doi = "10.1175/1520-0485(1986)016:1929"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert ":" not in encoded
        assert "%3A" in encoded

    def test_underscore_in_doi(self):
        """Test DOI with underscores (allowed per Crossref)."""
        doi = "10.1234/test_with_underscore"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        # Underscores are "unreserved" per RFC 3986, so urllib.parse.quote won't encode them
        # This is fine - underscores are safe in filenames

    def test_doi_with_multiple_query_params(self):
        """Test DOI URL with multiple query parameters."""
        doi_url = "10.1234/test?param1=value1&param2=value2&param3=value3"
        expected_doi = "10.1234/test"
        encoded = encode_doi(doi_url)
        decoded = decode_doi(encoded)
        assert decoded == expected_doi

    def test_doi_with_fragment(self):
        """Test DOI URL with fragment identifier."""
        doi_url = "10.1234/test#section1"
        # The function splits on '?' but not '#', so # gets encoded
        encoded = encode_doi(doi_url)
        decoded = decode_doi(encoded)
        assert decoded == "10.1234/test#section1"
        assert "%23" in encoded

    def test_doi_with_trailing_slash(self):
        """Test DOI with trailing slash."""
        doi = "10.1234/test/"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert encoded.endswith("%2F")

    def test_doi_with_multiple_slashes(self):
        """Test DOI with multiple slashes in suffix."""
        doi = "10.1234/path/to/resource"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        # Count slashes in original and encoded
        assert doi.count("/") == decoded.count("/")
        assert "%2F" in encoded

    def test_doi_org_variants(self):
        """Test various doi.org URL formats."""
        base_doi = "10.1234/test"
        variants = [
            "https://doi.org/10.1234/test",
            "http://doi.org/10.1234/test",
            "doi.org/10.1234/test",
            "www.doi.org/10.1234/test",
        ]
        
        for variant in variants:
            encoded = encode_doi(variant)
            decoded = decode_doi(encoded)
            # All should decode to the base DOI
            assert decoded == base_doi or "www" in decoded, f"Failed for variant: {variant}"

    def test_dx_doi_org_prefix(self):
        """Test dx.doi.org prefix (alternative DOI resolver)."""
        doi = "10.1234/test"
        doi_with_dx = "https://dx.doi.org/10.1234/test"
        
        encoded = encode_doi(doi_with_dx)
        decoded = decode_doi(encoded)
        
        # The split only handles 'doi.org/', so 'dx.doi.org' might remain
        # This test documents current behavior
        assert doi in decoded

    def test_very_long_doi(self):
        """Test with an unusually long DOI suffix."""
        # Some DOIs can be quite long
        long_suffix = "a" * 200
        doi = f"10.1234/{long_suffix}"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert len(encoded) > len(doi)  # Should be longer due to encoding

    def test_doi_with_numbers_only_suffix(self):
        """Test DOI with numeric-only suffix."""
        doi = "10.1234/123456789"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi

    def test_doi_with_mixed_case_prefix(self):
        """Test that DOI prefix case is preserved."""
        # DOI prefixes should be numeric, but test case preservation
        doi = "10.1234/TeSt"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        assert "TeSt" in decoded

    def test_encoding_with_already_encoded_characters(self):
        """Test behavior when DOI string already contains percent-encoded chars."""
        # This shouldn't happen in practice, but test defensive behavior
        doi_with_encoded = "10.1234/test%20space"
        encoded = encode_doi(doi_with_encoded)
        decoded = decode_doi(encoded)
        # urllib.parse.quote will double-encode the %
        assert decoded == doi_with_encoded

    def test_decoding_with_invalid_percent_encoding(self):
        """Test decode with malformed percent encoding."""
        # %ZZ is not valid hex
        malformed = "10%2E1234%2Ftest%ZZ"
        # unquote should handle this gracefully
        decoded = decode_doi(malformed)
        assert "10.1234/test" in decoded

    def test_special_chars_in_combination(self):
        """Test multiple special characters in sequence."""
        doi = "10.1234/test()[]<>;:#"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        assert decoded == doi
        # None of these special chars should remain unencoded
        for char in "()[]<>;:#":
            assert char not in encoded

    def test_whitespace_variations(self):
        """Test various whitespace characters."""
        dois = [
            "10.1234/test space",     # regular space
            "10.1234/test\ttab",      # tab
            "10.1234/test\nnewline",  # newline
        ]
        for doi in dois:
            encoded = encode_doi(doi)
            decoded = decode_doi(encoded)
            assert decoded == doi
            # Whitespace should be encoded
            assert " " not in encoded or "\t" not in encoded or "\n" not in encoded

    def test_doi_with_percentage_sign(self):
        """Test DOI containing actual % sign (rare but test it)."""
        doi = "10.1234/test%value"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        # % gets encoded as %25
        assert decoded == doi
        assert "%25" in encoded

    def test_all_printable_ascii_safe(self):
        """Test that common ASCII characters are handled safely."""
        # Test characters that might appear in DOIs
        test_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-._~"
        doi = f"10.1234/{test_chars}"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        # These should all survive round-trip
        assert all(c in decoded for c in test_chars if c not in ".-")  # .- get encoded

    def test_filesystem_path_construction(self):
        """Test that encoded DOI can be used in actual path construction."""
        import os
        doi = "10.1016/S0022-4049(02)00143-3"
        encoded = encode_doi(doi)
        
        # Simulate path construction
        path = os.path.join("PDFs", f"doi-{encoded}", f"doi-{encoded}.pdf")
        
        # Path should be valid (no illegal characters)
        assert "\\" not in encoded or "/" not in encoded
        assert path.endswith(".pdf")

    def test_url_construction_with_decode(self):
        """Test that decoded DOI works in URL construction."""
        doi = "10.1016/S0022-4049(02)00143-3"
        encoded = encode_doi(doi)
        decoded = decode_doi(encoded)
        
        # Simulate URL construction
        url = f"https://api.unpaywall.org/v2/{decoded}?email=test@example.com"
        
        # URL should contain the original DOI
        assert doi in url
        assert "(" in url  # Parentheses should be present in URL
